
extern crate serde_json;
extern crate serde;
extern crate rand;
extern crate chrono;

mod data;
mod s_state;
use data::*;
use s_state::*;
use std::{fs::{self, File}, env, collections::{HashMap, HashSet}};
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
    let data:String = fs::read_to_string(format!("../../problem.json/{}.json", index))?; 
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
    let mut best_score = -10000000000000000000.0;
    let mut best_name:String = String::new();

    let paths = fs::read_dir("../../solutions/").unwrap();
    for path in paths {
        let p = path.unwrap().path();
        let name = p.file_name().unwrap().to_str().unwrap();
        if !name.ends_with(".json") { continue; }
        if !name.starts_with(&format!("{}-", index)) { continue; }
        let mut submission = match File::open(format!("../../solutions/{}.submission.result", name))
        {
            Ok(file) => file,
            Err(err) => {
                //println!("error:{}", err);
                continue;
            }
        };
        let mut buf = String::new();
        submission.read_to_string(&mut buf)?;
        let mut score:f64 = match serde_json::from_str::<Submission>(&buf)
        {
            Ok(submission) => submission.Success.submission.score.Success,
            Err(err) => {
                0.0
            }
        };
        if name.contains("ueno") {
            score *= 1.006;
        }
        if score > best_score {
            best_name = name.to_string();
            best_score = score;
        }
    }
    let mut best_file = File::open(format!("../../solutions/{}", best_name))?;
    let mut buf = String::new();
    best_file.read_to_string(&mut buf)?;
    let mut best:Answer = serde_json::from_str(&buf)?;
    
    let mut swap_state = SwapState::new(&problem);
    let init_score = s_eval(&problem, &best.placements, &mut volumes, &mut swap_state);
    let musician_groups = swap_state.musician_groups;
    let active_attendees: HashSet<usize> = init_active_attendees(&problem, &best.placements);
    
    println!("{}:{} {}", index, init_score, best_score);
    let mut max_score = init_score;
    let mut max_result = best.placements;
    for j in 1..7
    {
        for i in 1..150 * j
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
                    pull_placements(&mut placements, &problem, 0.45 * rng.gen_range(0.1..1.0), &musician_groups, &active_attendees, &mut rng);
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
                println!("up: {} {} {},{}:{}:{} {}", index, best_name, j, i, pulls, score, (score - init_score) / init_score);
            }
        }
    }
    if max_score > init_score * 1.0035
    {
        let score = s_eval(&problem, &max_result, &mut volumes, &mut &mut SwapState::new(&problem));
        println!("end: {} {}:{} {} ", index, best_name, score, (score - init_score) / init_score);
        let answer:Answer = Answer { placements:max_result, volumes };
        let answer_string = serde_json::to_string(&answer)?;
        while best_name.len() > 30
        {
            best_name.remove(0); 
        }
        let name = format!("shohei15-6-{}-{}", seed, best_name);
        fs::write(
            format!("../../solutions/{}-{}", index, name), 
            &answer_string
        )?;
        fs::write(
            format!("../../solutions/{}-{}.myeval", index, name), 
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
fn pull_placements<R:Rng>(
    placements:&mut Vec<Point>,
    problem:&Problem,
    power:f64,
    musician_groups:&HashMap<usize, Vec<usize>>,
    active_attendees:&HashSet<usize>,
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

            for i in active_attendees.iter()
            {
                let attendee = &problem.attendees[*i];
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
                    vx += dx / d * v / 3.0;
                    vy += dy / d * v / 3.0;
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
