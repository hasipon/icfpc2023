
extern crate serde_json;
extern crate serde;
extern crate rand;
extern crate chrono;

mod data;
mod s_state;
use data::*;
use s_state::*;
use std::{fs::{self, File}, env};
use rand::{Rng, rngs::ThreadRng};
use chrono::prelude::*;
use std::io::Read;


fn main() -> Result<(), Box<dyn std::error::Error>> {
    let timestamp = Utc::now().timestamp();

    let args: Vec<String> = env::args().collect();
    let id = if args.len() <= 1 { "9" } else { &args[1] };
    solve(id, timestamp)?;

    Ok(())
}

fn solve(index:&str, timestamp:i64) -> Result<(), Box<dyn std::error::Error>> {
    let data:String = fs::read_to_string(format!("../../../problem.json/{}.json", index))?;
    let mut problem:Problem = serde_json::from_str(&data)?;
    let x = problem.stage_bottom_left.0 + 10.0;  
    let y = problem.stage_bottom_left.1 + 10.0;
    let w = problem.stage_width - 20.0;
    let h = problem.stage_height - 20.0;
    let mut seed = rand::thread_rng().gen_range(1..10000);
    let mut rng: rand::rngs::StdRng = rand::SeedableRng::seed_from_u64(seed);
    
    let mut volumes = vec![1.0; problem.musicians.len()];
    let iindex:i64 = index.parse()?;
    if iindex > 55 {
        problem.extention = Option::Some(());
    }

    // 最善解をとってくる
    let mut best_file = File::open(format!("../../../solutions/{}-shiota13.json", index))?;
    let mut buf = String::new();
    best_file.read_to_string(&mut buf)?;
    let mut best:Answer = serde_json::from_str(&buf)?;
    let init_score = s_eval(&problem, &best.placements, &mut volumes);
    let mut max_score = init_score;
    let mut max_result = best.placements;

    while(true)
    {
        let mut found = false;
        for i in 1..180
        {
            let mut placements = max_result.clone();
            let swaps = rng.gen_bool(0.04);
            if !swaps || rng.gen_bool(0.9) 
            {
                while !randomize(&mut placements, &problem, 600.0 / (i * 2) as f64, (i * 2) as f64, &mut rng)
                {
                }
                if rng.gen_bool(0.02) {
                    pull_placements(&mut placements, &problem, 0.1);
                }
                for _ in 0..2000
                {
                    if separate_placements(x, y, w, h, &mut placements, &mut rng) {
                        break;
                    }
                }
            }
            if swaps {
                try_swap(&problem, &mut placements, &mut rng, &mut SwapState::new(&problem));
            }
            
            let score = s_eval(&problem, &placements, &mut volumes);
            if score > max_score
            {
                max_score = score;
                max_result = placements;
                eprintln!("max_score: {}", max_score);
                found = true
            }
        }
        if !found {
            break;
        }
    }
    if max_score > init_score * 1.001
    {
        let score = s_eval(&problem, &max_result, &mut volumes);
        let answer:Answer = Answer { placements:max_result, volumes };
        let answer_string = serde_json::to_string(&answer)?;
        let name = format!("shohei15-2-{}-{}", seed, "shiota-13");
        fs::write(
            format!("../../../solutions/{}-{}.json", index, name),
            &answer_string
        )?;
        fs::write(
            format!("../../../solutions/{}-{}.myeval", index, name),
            &score.to_string()
        )?;
    }
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
fn pull_placements(
    placements:&mut Vec<Point>,
    problem:&Problem,
    power:f64) {
    let mut max_v = 0.0;
    let mut vs = Vec::new();
    for i in 0..placements.len()
    {
        let mut p1 = placements[i];
        let index = problem.musicians[i] as usize;
        let mut vx = 0.0;
        let mut vy = 0.0;
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
        if v > max_v {
            max_v = v;
        }
        vs.push((vx, vy));
    }
    if max_v == 0.0 { return; }
    for (i, v) in vs.iter().enumerate()
    {
        placements[i].x += v.0 / max_v * power;
        placements[i].y += v.1 / max_v * power;
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
