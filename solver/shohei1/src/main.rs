
extern crate serde_json;
extern crate serde;
extern crate rand;
extern crate chrono;

mod data;
use data::*;
use std::fs;
use rand::Rng;
use chrono::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let timestamp = Utc::now().timestamp();
    for i in 31..45 {
        solve(i + 1, timestamp)?;
    }
    Ok(())
}

fn solve(index:u32, timestamp:i64) -> Result<(), Box<dyn std::error::Error>> {
    let data:String = fs::read_to_string(format!("../../problem.json/{}.json", index))?; 
    let problem:Problem = serde_json::from_str(&data)?;
    let mut placements = Vec::new();
    let x = problem.stage_bottom_left.0 + 10.0;  
    let y = problem.stage_bottom_left.1 + 10.0;
    let w = problem.stage_width - 20.0;
    let h = problem.stage_height - 20.0;
    let mut max_taste = 0.0;
    let mut rng = rand::thread_rng();
    for _ in &problem.musicians
    {
        placements.push(
            Point{
                x: x + rng.gen_range(0.0..w + 0.00001), 
                y: y + rng.gen_range(0.0..h + 0.00001)
            }
        );
    }
    for attendee in &problem.attendees
    {
        for taste in &attendee.tastes
        {
            if max_taste < taste.abs()
            {
                max_taste = taste.abs();
            }
        }
    }
    for i in 0..100
    {
        let power = 30.0 / ((i as f64).powf(1.4)  + 2.0) / max_taste; 
        pull_placements(&mut placements, &problem, power);
        separate_placements(x, y, w, h, &mut placements, &mut rng);
    }
    for _ in 0..2000
    {
        if separate_placements(x, y, w, h, &mut placements, &mut rng) {
            break;
        }
    }

    let answer:Answer = Answer { placements };
    let answer_string = serde_json::to_string(&answer)?;
    fs::write(
        format!("../../solutions/{}-shohei1-{}.json", index, timestamp), 
        &answer_string
    )?;
    Ok(())
}

fn pull_placements(
    placements:&mut Vec<Point>,
    problem:&Problem,
    power:f64) {
    
    for i in 0..placements.len()
    {
        let mut p1 = placements[i];
        let index = problem.musicians[i] as usize;
        for attendee in &problem.attendees
        {
            let dx = p1.x - attendee.x;
            let dy = p1.y - attendee.y;
            let d2 = dx * dx + dy * dy;
            let value = power * attendee.tastes[index] / d2.sqrt();
            p1.x -= dx * value;
            p1.y -= dy * value;
        }
        placements[i] = p1;
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
            p1.x = x + rng.gen_range(0.0..0.0000001);
            hit = true;
        }
        if p1.y < y
        {
            p1.y = y + rng.gen_range(0.0..0.0000001);
            hit = true;
        }
        if p1.x > x + w
        {
            p1.x = x + w - rng.gen_range(0.0..0.0000001);
            hit = true;
        }
        if p1.y > y + h
        {
            p1.y = y + h - rng.gen_range(0.0..0.0000001);
            hit = true;
        }
        for j in (i + 1)..placements.len()
        {
            let mut p2 = placements[j];
            
            let dx = p1.x - p2.x;
            let dy = p1.y - p2.y;
            let d2 = dx * dx + dy * dy;
            if d2 < 400.0
            {
                let d = d2.sqrt();
                let power = (20.001 - d) * 0.8;
                p1.x += power * (dx / d);
                p2.x -= power * (dx / d);
                p1.y += power * (dy / d);
                p2.y -= power * (dy / d);
                hit = true;
            }
            placements[j] = p2;
        }
        if p1.x < x
        {
            p1.x = x;
            hit = true;
        }
        if p1.y < y
        {
            p1.y = y;
            hit = true;
        }
        if p1.x > x + w
        {
            p1.x = x + w;
            hit = true;
        }
        if p1.y > y + h
        {
            p1.y = y + h;
            hit = true;
        }
        placements[i] = p1;
    }
    !hit
}
