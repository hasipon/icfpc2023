
extern crate serde_json;
extern crate serde;
extern crate rand;
extern crate chrono;
extern crate core;

mod data;
mod state;

use data::*;
use state::*;
use std::{fs, cmp::Ordering, env, process};
use rand::Rng;
use chrono::prelude::*;
use std::collections::BinaryHeap;
use std::process::exit;

use crate::state::CacheState;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let timestamp = Utc::now().timestamp();

    let args: Vec<String> = env::args().collect();
    let id = match env::var("PROBLEM_ID") {
        Ok(val) => val,
        Err(err) => {
            println!("{}: {}", err, "PROBLEM_ID");
            process::exit(1);
        },
    };
    solve(&id, timestamp)?;
    Ok(())
}

fn calcPerm(problem:&Problem) -> Vec<usize> {
    let mut forSort  = Vec::new();
    for(i, t) in problem.musicians.iter().enumerate() {
        forSort.push((t, i));
    }
    forSort.sort();
    let mut perm = Vec::new();
    for v in forSort {
        perm.push(v.1);
    }
    perm
}

fn coreSolveGaishu(problem:&Problem, offsetX:f64, offsetY:f64) -> (f64, Vec<Point>, Vec<f64>) {
    let mut cache = CacheState::new(&problem);
    let mut rng: rand::rngs::StdRng = rand::SeedableRng::seed_from_u64(3);
    let mut placements = Vec::new();
    let mut volume = Vec::new();
    let x = problem.stage_bottom_left.0 + 10.0 + offsetX;
    let y = problem.stage_bottom_left.1 + 10.0 + offsetY;
    let w = problem.stage_width - 20.0;
    let h = problem.stage_height - 20.0;

    // 外側に配置する
    {
        let mut dx = 0.0;
        let mut dy: f64 = -1.0;
        let mut left = x;
        let mut right = x + w;
        let mut top = y;
        let mut bottom = y + h;
        let mut cx = left;
        let mut cy = bottom;

        let perm = calcPerm(&problem);
        let mut permI = 0;
        while permI < problem.musicians.len()
        {
            if (bottom - top).abs() < 10. {
                return (0.0, Vec::new(), Vec::new())
            }
            if (left - right).abs() < 10. {
                return (0.0, Vec::new(), Vec::new())
            }
            let i = perm[permI];
            permI += 1;
            let speed = 10.0;
            placements.push(Point{x:cx, y:cy});
            volume.push(0.);

            cx += dx * speed;
            cy += dy * speed;
            if cy > bottom {
                right -= 10.0;
                cy = bottom;
                cx = right;

                dx = -1.0;
                dy = 0.0;
            }
            if cx < left {
                bottom -= 10.0;
                cx = left;
                cy = bottom;

                dy = -1.0;
                dx = 0.0;
            }
            if cy < top {
                left += 10.0;
                cy = top;
                cx = left;

                dx = 1.0;
                dy = 0.0;
            }
            if cx > right {
                top += 10.0;
                cx = right;
                cy = top;

                dx = 0.0;
                dy = 1.0;
            }
        }
    }

    let mut lastI  = 0.;
    for i in 0..20
    {
        let now = eval(&problem, &placements, &mut volume, &mut cache);
        if now == lastI{
            break;
        }
        lastI = now;
        try_bulk_swap(&problem, &mut placements, &mut volume,&mut rng, &mut cache);
        let mut lastJ  = 0.;
        for j in 0..10
        {
            let now = eval(&problem, &placements, &mut volume, &mut cache);
            if now == lastJ{
                break;
            }
            lastJ = now;
            try_swap(&problem, &mut placements, &mut volume, &mut rng, &mut cache);
        }
    }
    let score = eval(&problem, &placements, &mut volume,&mut cache);
    return (score, placements, volume);
}

fn coreSolveHeat(problem:&Problem, heat:&Vec<(f64,f64)>, index:i32, offsetX:f64, offsetY:f64) -> (f64, Vec<Point>, Vec<f64>) {
    let mut cache = CacheState::new(&problem);
    let mut rng: rand::rngs::StdRng = rand::SeedableRng::seed_from_u64(3);
    let mut placements = Vec::new();
    let mut volume = Vec::new();

    let xRange = (problem.stage_bottom_left.0 + 10.0 , problem.stage_bottom_left.0 + problem.stage_width - 10.0);
    let yRange = (problem.stage_bottom_left.1 + 10.0 , problem.stage_bottom_left.1 + problem.stage_height - 10.0);
    let mut amari = Vec::new();
    for (xx, yy) in heat {
        let x = offsetX + xx;
        let y = offsetY + yy;
        let hoge = 5.;
        if (x - xRange.0).abs() > hoge && (x - xRange.1).abs() > hoge && (y - yRange.0).abs() > hoge && (y - yRange.1).abs() >hoge {
            amari.push((x,y));
            continue;
        }
        placements.push(Point{x, y});
        volume.push(0.);
        if placements.len() >= problem.musicians.len() {
            break;
        }
    }
    while placements.len() != problem.musicians.len()
    {
        let (x, y) = amari.pop().unwrap();
        placements.push(Point{x, y});
        volume.push(0.);
    }

    let mut lastI  = 0.;
    while(true)
    {
        let now = eval(&problem, &placements, &mut volume, &mut cache);
        if now == lastI{
            break;
        }
        lastI = now;
        try_bulk_swap(&problem, &mut placements, &mut volume, &mut rng, &mut cache);
        let mut lastJ  = 0.;
        for j in 0..10
        {
            let now = eval(&problem, &placements, &mut volume, &mut cache);
            if now == lastJ{
                break;
            }
            lastJ = now;
            try_swap(&problem, &mut placements, &mut volume, &mut rng, &mut cache);
        }
    }
    let score = eval(&problem, &placements, &mut volume, &mut cache);
    eprintln!("{}:{}", index, score);
    return (score, placements, volume);
}

fn solve(index:&str, timestamp:i64) -> Result<(), Box<dyn std::error::Error>> {
    let repoRoot = match env::var("REPO_ROOT") {
        Ok(val) => val,
        Err(err) => {
            println!("{}: {}", err, "REPO_ROOT");
            process::exit(1);
        },
    };
    let data:String = fs::read_to_string(format!("{}/problem.json/{}.json", repoRoot, index))?;
    let mut problem:Problem = serde_json::from_str(&data)?;

    let heatStr:String = fs::read_to_string(format!("{}/solver/ueno4/{}.csv", repoRoot, index))?;
    let mut heat = Vec::new();
    for (a, line) in heatStr.split("\n").enumerate(){
        let sp:Vec<&str> = line.split(',').collect();
        if sp.len () < 3 {
            break;
        }
        let x = sp[1].parse::<f64>()?;
        let y = sp[2].parse::<f64>()?;
        heat.push((x, y));
    }

    let index:i32 = index.parse()?;
    if index > 55 {
        problem.extention = Option::Some(());
    }
    let mut bestScore = 0.0;
    let mut bestPlacements = Vec::new();
    let mut bestVolumes = Vec::new();
    for offsetX in [10 - problem.stage_bottom_left.0 as i64 % 10, 10 - (problem.stage_bottom_left.0 + problem.stage_width) as i64 % 10]
    {
        for offsetY in  [10- problem.stage_bottom_left.1 as i64 % 10, 10 - (problem.stage_bottom_left.1 + problem.stage_height) as i64 % 10]{
            let (score, placements, volumes) = coreSolveHeat(&problem, &heat, index, (offsetX % 10) as f64, (offsetY %10) as f64);
            eprintln!("heat {}", score);
            if score > bestScore {
                bestScore = score;
                bestPlacements = placements;
                bestVolumes = volumes;
            }
        }
    }
    for offsetX in [0, (10 + (problem.stage_bottom_left.0) as i64 %10 - (problem.stage_bottom_left.0 + problem.stage_width) as i64 % 10) %10] {
        for offsetY in [0, (10 + (problem.stage_bottom_left.1) as i64 %10 - (problem.stage_bottom_left.1 + problem.stage_height) as i64 % 10) %10] {
            let (score, placements, volumes) = coreSolveGaishu(&problem, offsetX as f64, offsetY as f64);
            eprintln!("gaishu {}", score);
            if score > bestScore {
                bestScore = score;
                bestPlacements = placements;
                bestVolumes = volumes;
            }
        }
    }
    eprintln!("{}\t{}", bestScore, index);

    let answer:Answer = Answer { placements: bestPlacements, volumes: bestVolumes};
    let answer_string = serde_json::to_string(&answer)?;
    let name = "shiota10";
    println!("{}", &answer_string);
    // fs::write(
    //     format!("../../solutions/{}-{}.json", index, name),
    //     &answer_string
    // )?;
    // fs::write(
    //     format!("../../solutions/{}-{}-{}.myeval", index, name, timestamp),
    //     &bestScore.to_string()
    // )?;

    Ok(())
}

fn try_bulk_swap<R:Rng>(problem:&Problem, placements:&mut Vec<Point>, volume:&mut Vec<f64>, rng:&mut R, cache:&mut CacheState) {
    let tasteN = problem.attendees[0].tastes.len();
    for ti in 0 .. tasteN {
        for tj in ti + 1 .. tasteN {
            let mut scoreBefore = 0.;
            let mut scoreAfter= 0.;
            let mut swapHist = Vec::new();
            for i in 0 .. placements.len() {
                if problem.musicians[i] != ti {
                    continue;
                }
                let mut found = false;
                for j in i + 1 .. placements.len(){
                    if problem.musicians[j] != tj {
                        continue;
                    }
                    scoreBefore = scoreBefore + eval_placement(problem, placements, volume, i, false, cache) + eval_placement(problem, placements, volume, j, false, cache);
                    let scorei = cache.placement_score[i];
                    let scorej = cache.placement_score[j];
                    swapHist.push((i, j, scorei, scorej));

                    placements.swap(i, j);
                    cache.swap_placement(i, j);
                    found = true;
                    scoreAfter =  scoreAfter +  eval_placement(problem, placements, volume,i, false, cache) + eval_placement(problem, placements, volume, j, false, cache);
                    break;
                }
                if found {
                    continue;
                }
                for j in i + 1 .. placements.len(){
                    if problem.musicians[j] == ti {
                        continue;
                    }
                    scoreBefore = scoreBefore + eval_placement(problem, placements, volume, i, false, cache) + eval_placement(problem, placements, volume, j, false, cache);
                    let scorei = cache.placement_score[i];
                    let scorej = cache.placement_score[j];
                    swapHist.push((i, j, scorei, scorej));

                    placements.swap(i, j);
                    cache.swap_placement(i, j);
                    found = true;
                    scoreAfter =  scoreAfter +  eval_placement(problem, placements, volume, i, false, cache) + eval_placement(problem, placements, volume, j, false, cache);
                    break;
                }
            }
            if scoreBefore > scoreAfter {
                swapHist.reverse();
                for (i, j, scorei, scorej) in swapHist {
                    placements.swap(i, j);
                    cache.swap_placement(i, j);
                    cache.placement_score[i] = scorei;
                    cache.placement_score[j] = scorej;
                }
            }
        }
    }
}

fn try_swap<R:Rng>(problem:&Problem, placements:&mut Vec<Point>, volume:&mut Vec<f64>, rng:&mut R, cache:&mut CacheState) {
    for i in 0..placements.len()
    {
        for j in i + 1..placements.len()
        {
            if problem.musicians[i] == problem.musicians[j] {
                continue;
            }
                let score1 = eval_placement(problem, placements,volume, i, false, cache) + eval_placement(problem, placements, volume, j, false, cache);
                let scorei = cache.placement_score[i];
                let scorej = cache.placement_score[j];

                placements.swap(i, j);
                cache.swap_placement(i, j);
                
                let score2 = eval_placement(problem, placements, volume,i, false, cache) + eval_placement(problem, placements, volume, j, false, cache);
                if score1 > score2 {
                    placements.swap(i, j);
                    cache.swap_placement(i, j);
                    cache.placement_score[i] = scorei;
                    cache.placement_score[j] = scorej;
            }
        }
    } 
}

fn eval(problem:&Problem, placements:&Vec<Point>, volume:&mut Vec<f64>, cache:&mut CacheState) -> f64 {
    let mut result = 0.0;
    for i in 0..placements.len()
    {
        result += eval_placement(problem, placements, volume, i, true, cache);
    }
    result
}

// 各ミュージシャンごとの観客の評価値の合算
fn eval_placement(problem:&Problem, placements:&Vec<Point>, volume:&mut Vec<f64>, index:usize, total:bool, cache:&mut CacheState) -> f64 {
    let mut nearest_d = 10000000000000000000000.0;
    let mut nearest_dir = 0.0;
    let center = placements[index];

    if cache.placement_sights[index].is_none()
    {
        let mut nodes = Vec::new();
        // ミュージシャンのふちを追加
        for (i, p) in placements.iter().enumerate()
        {
            if *p == center || i == index { continue; }
            let dx = p.x - center.x;
            let dy = p.y - center.y;
            let d = (dx * dx + dy * dy).sqrt();
            let asin = if d <= 5.0 { 
                std::f64::consts::PI / 2.0 
            } else {
                (5.0 / d).asin()
            };
            let dir = dx.atan2(dy) - asin; 
            nodes.push(
                EvalNode {
                    dir,
                    d,
                    kind:EvalNodeKind::Start(asin * 2.0),
                }
            );
            if d < nearest_d {
                nearest_d = d;
                nearest_dir = dir;
            }
        }
        // 柱のふちを追加
        for p in &problem.pillars
        {
            let dx = p.center.0 - center.x;
            let dy = p.center.1 - center.y;
            let d = (dx * dx + dy * dy).sqrt();
            let asin = if d <= p.radius { 
                std::f64::consts::PI / 2.0 
            } else {
                (p.radius / d).asin()
            };
            let dir = dx.atan2(dy) - asin; 
            nodes.push(
                EvalNode {
                    dir,
                    d,
                    kind:EvalNodeKind::Start(asin * 2.0),
                }
            );
            if d < nearest_d {
                nearest_d = d;
                nearest_dir = dir;
            }
        }
        // 客の中心を追加
        for (i, a) in problem.attendees.iter().enumerate()
        {
            let dx = a.x - center.x;
            let dy = a.y - center.y;
            let d = (dx * dx + dy * dy).sqrt();
            let dir = dx.atan2(dy); 
            nodes.push(
                EvalNode {
                    dir,
                    d,
                    kind:EvalNodeKind::Center(i),
                }
            );
        }
        for node in &mut nodes
        {
            node.dir -= nearest_dir;
            if node.dir < 0.0
            {
                node.dir += std::f64::consts::PI * 2.0;
            }
        }
        nodes.sort_unstable_by(|a, b| (a.dir, a.d).partial_cmp(&(b.dir, b.d)).unwrap_or(Ordering::Equal));    
        
        let mut heap:BinaryHeap<EvalNode> = BinaryHeap::new();
        let mut sights = Vec::new();
        for node in &nodes
        {
            let next_dir = node.dir;
            while !heap.peek().is_none() && heap.peek().unwrap().end_dir() < next_dir
            {
                heap.pop();
            }
            match &node.kind
            {
                &EvalNodeKind::Start(_) => { heap.push(node.clone()); }
                &EvalNodeKind::Center(aindex) => {
                    if heap.peek().is_none() || node.d < heap.peek().unwrap().d {
                        sights.push(Sight{attendee:aindex, d2: node.d * node.d });
                    }
                },
            }
        }
        cache.placement_sights[index] = Option::Some(sights);
    }

    let mut result = 0.0;
    for sight in cache.placement_sights[index].as_ref().unwrap() {
        let a = &problem.attendees[sight.attendee];
        result += 1000000.0 * a.tastes[problem.musicians[index]] / sight.d2;
    }
    cache.placement_score[index] = result;

    if problem.extention.is_some() {
        let mut q = 1.0;
        let mut bonus = 0.0;
        for i in cache.musician_groups.get(&problem.musicians[index]).unwrap()
        {
            if *i == index { continue; }
            let p = placements[*i];
            let dx = center.x - p.x;
            let dy = center.y - p.y;
            let d = (dx * dx + dy * dy).sqrt();
            let value = 1.0 / d.max(10.0);
            q += value;
            bonus += value * cache.placement_score[*i];
        }
        result *= q;
        
        if !total
        {
            result += bonus;
        }
    }
    if result > 0. {
        volume[index] = 10.;
    } else {
        volume[index] = 0.;
    }
    result * volume[index]
}

