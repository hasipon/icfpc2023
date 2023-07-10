
extern crate serde_json;
extern crate serde;
extern crate rand;
extern crate chrono;

mod data;
mod s_state;
mod g_state;
use data::*;
use s_state::*;
use g_state::*;
use std::{fs, env, collections::HashMap};
use rand::Rng;
use chrono::prelude::*;


fn main() -> Result<(), Box<dyn std::error::Error>> {
    let timestamp = Utc::now().timestamp();

    let args: Vec<String> = env::args().collect();
    let id = if args.len() <= 1 { "2" } else { &args[1] };
    solve(id, timestamp)?;

    Ok(())
}

fn solve(index:&str, timestamp:i64) -> Result<(), Box<dyn std::error::Error>> {
    let data:String = fs::read_to_string(format!("../../problem.json/{}.json", index))?; 
    let mut problem:Problem = serde_json::from_str(&data)?;
    let x = problem.stage_bottom_left.0 + 10.0;  
    let y = problem.stage_bottom_left.1 + 10.0;
    let w = problem.stage_width - 20.0;
    let h = problem.stage_height - 20.0;
    let seed = rand::thread_rng().gen_range(0..20000);
    let mut rng: rand::rngs::StdRng = rand::SeedableRng::seed_from_u64(seed);
    
    let mut volumes = vec![1.0; problem.musicians.len()];
    let iindex:i64 = index.parse()?;
    if iindex > 55 {
        problem.extention = Option::Some(());
    }

    // ランダムに配置する
    let mut grid_state = GridState::new(&problem);
    let mut placements = grid_state.init_grid(&problem, &mut rng);
    for center in &placements
    {
        if center.x.is_nan() { panic!("x NaN!!!!!! {} {:?}", index, placements); }
        if center.y.is_nan() { panic!("y NaN!!!!!! {} {:?}", index, placements); }
    }

    println!("{}", placements.len());
    let mut swap_state = SwapState::new(&problem);
    for i in 0..1
    {
        println!("{}: s{}:{}", index, i, s_eval(&problem, &placements, &mut volumes, &mut swap_state));
        try_swap(&problem, &mut placements, &mut rng, &mut swap_state);
        swap_state = SwapState::new(&problem);
    
        println!("{}: g{}:{}", index, i, s_eval(&problem, &placements, &mut volumes, &mut swap_state));
    }
    
    let init_score = s_eval(&problem, &placements, &mut volumes, &mut SwapState::new(&problem));
    let mut max_score = init_score;
    let mut max_result = placements;
    let musician_groups = SwapState::new(&problem).musician_groups;

    for j in 1..6
    {
        for i in 1..125 * j
        {
            let mut placements = max_result.clone();
            let swaps = rng.gen_bool(0.01);
            let mut pulls = false;
            if !swaps || rng.gen_bool(0.8) 
            {
                pulls = rng.gen_bool(0.5);
                if !pulls || rng.gen_bool(0.5) 
                {
                    while !randomize(&mut placements, &problem, 600.0 / (i * 2) as f64, (i * 2) as f64, &mut rng)
                    {
                    }
                }
                
                if pulls {
                    pull_placements(&mut placements, &problem, 0.45 * rng.gen_range(0.1..1.0), &musician_groups, &mut rng);
                }
                for _ in 0..2000
                {
                    if separate_placements(x, y, w, h, &mut placements, &mut rng) {
                        break;
                    }
                }
            }
            let mut swap_state = SwapState::new(&problem);
            if swaps {
                try_swap(&problem, &mut placements, &mut rng, &mut swap_state);
            }
            
            let score = s_eval(&problem, &placements, &mut volumes, &mut swap_state);
            if score > max_score
            {
                max_score = score;
                max_result = placements;
                println!("up: {} {},{}:{}:{} {}", index, j, i, pulls, score, (score - init_score) / init_score);
            }
        }
    }

    let placements = max_result;
    let score = s_eval(&problem, &placements, &mut volumes, &mut swap_state);
    println!("end:{}:{}", index, score);

    let answer:Answer = Answer { placements, volumes };
    let answer_string = serde_json::to_string(&answer)?;
    let name = format!("shohei17-{}", seed);
    fs::write(
        format!("../../solutions/{}-{}.json", index, name), 
        &answer_string
    )?;
    fs::write(
        format!("../../solutions/{}-{}.myeval", index, name), 
        &score.to_string()
    )?;
    Ok(())
}

fn randomize<R:Rng>(
    placements:&mut Vec<Point>,
    problem:&Problem,
    power:f64,
    rate:f64,
    rng:&mut R) -> bool {

    let mut changed = false;
    let rate = (rng.gen_range(0.000001..0.0003) * rate).min(1.0);
    for i in 0..placements.len()
    {
        if rng.gen_bool(rate) {
            placements[i].x += rng.gen_range(-power..power);
            placements[i].y += rng.gen_range(-power..power);
            changed = true;
        }
    }
    changed
}
fn pull_placements<R:Rng>(
    placements:&mut Vec<Point>,
    problem:&Problem,
    power:f64,
    musician_groups:&HashMap<usize, Vec<usize>>,
    rng:&mut R) {

    let mut max_v = 0.0;
    let mut vs = Vec::new();
    let grouping = problem.extention.is_some() && rng.gen_bool(0.5);
    for i in 0..placements.len()
    {
        let mut p1 = placements[i];
        let index = problem.musicians[i] as usize;
        let mut vx = 0.0;
        let mut vy = 0.0;
        if rng.gen_bool(0.03) {

            for attendee in &problem.attendees
            {
                let dx = p1.x - attendee.x;
                let dy = p1.y - attendee.y;
                let d2 = dx * dx + dy * dy;
                let speed = attendee.tastes[index] / d2;
                vx -= dx * speed;
                vy -= dy * speed;
            }

            let v = (vx * vx + vy * vy).sqrt();

            if grouping {
                for j in musician_groups.get(&problem.musicians[i]).unwrap()
                {
                    if *j == i { continue; }
                    if rng.gen_bool(0.5) { continue; }
                    let p2 = placements[*j];
                    let dx = p2.x - p1.x;
                    let dy = p2.y - p1.y;
                    let d = (dx * dx + dy * dy).sqrt();
                    vx += dx / d * v / 2.0;
                    vy += dy / d * v / 2.0;
                }
            }
        }
        let rate = rng.gen_range(0.0..1.0);
        vx *= rate;
        vy *= rate;
        let v = (vx * vx + vy * vy).sqrt();
        if v > max_v {
            max_v = v;
        }
        vs.push((vx, vy));
    }
    if max_v == 0.0 { return; }
    for (i, v) in vs.iter().enumerate()
    {
        let rate = rng.gen_range(0.0..1.0);
        placements[i].x += v.0 / max_v * power * rate;
        placements[i].y += v.1 / max_v * power * rate;
    }
}


fn separate_placements<R:Rng>(
    x:f64,
    y:f64,
    w:f64,
    h:f64,
    placements:&mut Vec<Point>, 
    rng:&mut R) -> bool {
    let mut hit = false;
    for i in 0..placements.len()
    {
        let mut p1 = placements[i];
        if p1.x < x
        {
            let d = p1.x - x;
            p1.x = x + d / 10.0;
            hit = true;
        }
        if p1.y < y
        {
            let d = p1.y - y;
            p1.y = y + d / 10.0;
            hit = true;
        }
        if p1.x > x + w
        {
            let d = p1.x - (x + w);
            p1.x = x + w + d / 10.0;
            hit = true;
        }
        if p1.y > y + h
        {
            let d = p1.y - (y + h);
            p1.y = y + h + d / 10.0;
            hit = true;
        }
        for j in (i + 1)..placements.len()
        {
            let mut p2 = placements[j];
            
            let dx = p1.x - p2.x;
            let dy = p1.y - p2.y;
            let d2 = dx * dx + dy * dy;
            if d2 < 100.0
            {
                let d = d2.sqrt() + 0.0001;
                let power = (10.2 - d) * 0.48;
                p1.x += power * (dx / d);
                p2.x -= power * (dx / d);
                p1.y += power * (dy / d);
                p2.y -= power * (dy / d);
                hit = true;
            }
            placements[j] = p2;
        }
        //if p1.x < x
        //{
        //    p1.x = x;
        //    hit = true;
        //}
        //if p1.y < y
        //{
        //    p1.y = y;
        //    hit = true;
        //}
        //if p1.x > x + w
        //{
        //    p1.x = x + w;
        //    hit = true;
        //}
        //if p1.y > y + h
        //{
        //    p1.y = y + h;
        //    hit = true;
        //}
        placements[i] = p1;
    }
    !hit
}
