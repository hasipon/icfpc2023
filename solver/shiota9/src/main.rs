
extern crate serde_json;
extern crate serde;
extern crate rand;
extern crate chrono;
extern crate core;

mod data;
mod state;

use data::*;
use state::*;
use std::{fs, cmp::Ordering, env};
use rand::Rng;
use chrono::prelude::*;
use std::collections::BinaryHeap;

use crate::state::CacheState;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let timestamp = Utc::now().timestamp();

    let args: Vec<String> = env::args().collect();
    let id = if args.len() <= 1 { "81" } else { &args[1] };
    solve(id, timestamp)?;
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

fn coreSolve(problem:&Problem, w:f64, h:f64, index:i32) -> (f64, Vec<Point>) {
    let mut cache = CacheState::new(&problem);
    let mut rng: rand::rngs::StdRng = rand::SeedableRng::seed_from_u64(3);
    let mut placements = Vec::new();
    let x = problem.stage_bottom_left.0 + 10.0;
    let y = problem.stage_bottom_left.1 + 10.0;
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
                return (0.0, Vec::new())
            }
            if (left - right).abs() < 10. {
                return (0.0, Vec::new())
            }
            let i = perm[permI];
            permI += 1;
            let speed = 10.0;
            placements.push(Point{x:cx, y:cy});

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
        let now = eval(&problem, &placements, &mut cache);
        if now == lastI{
            break;
        }
        lastI = now;
        try_bulk_swap(&problem, &mut placements, &mut rng, &mut cache);
        let mut lastJ  = 0.;
        for j in 0..10
        {
            let now = eval(&problem, &placements, &mut cache);
            if now == lastJ{
                break;
            }
            lastJ = now;
            try_swap(&problem, &mut placements, &mut rng, &mut cache);
        }
    }
    let score = eval(&problem, &placements, &mut cache);
    return (score, placements);
}

fn solve(index:&str, timestamp:i64) -> Result<(), Box<dyn std::error::Error>> {
    let data:String = fs::read_to_string(format!("../../problem.json/{}.json", index))?; 
    let mut problem:Problem = serde_json::from_str(&data)?;
    let index:i32 = index.parse()?;
    if index > 55 {
        problem.extention = Option::Some(());
    }
    let w = problem.stage_width - 20.0;
    let h = problem.stage_height - 20.0;

    let mut bestScore = 0.0;
    let mut bestPlacements = Vec::new();
    let unit = 10.;
    let mut nowW = 0.0;
    while nowW < w {
        let mut nowH = 0.;
        while nowH < h {
            let (score, placements) = coreSolve(&problem, nowW, nowH, index);
            if score > bestScore {
                bestScore = score;
                bestPlacements = placements;
                println!("{}\t{}\t{}", bestScore, nowW, nowH);
            }
            nowH = nowH + unit;
        }
        nowW = nowW + unit;
    }
    println!("{}\t{}", bestScore, index);

    let answer:Answer = Answer { placements: bestPlacements };
    let answer_string = serde_json::to_string(&answer)?;
    let name = "shiota9";
    fs::write(
        format!("../../solutions/{}-{}-{}.json", index, name,timestamp),
        &answer_string
    )?;
    fs::write(
        format!("../../solutions/{}-{}-{}.myeval", index, name, timestamp),
        &bestScore.to_string()
    )?;

    Ok(())
}

fn try_bulk_swap<R:Rng>(problem:&Problem, placements:&mut Vec<Point>, rng:&mut R, cache:&mut CacheState) {
    let tasteN = problem.attendees[0].tastes.len();
    for ti in 0 .. tasteN {
        for tj in ti + 1 .. tasteN {
            let mut scoreBefore = 0.;
            let mut scoreAfter= 0.;
            let mut swapHist = Vec::new();
            let beforeScoreStrict = eval(problem, placements, cache);
            for i in 0 .. placements.len() {
                if problem.musicians[i] != ti {
                    continue;
                }
                let mut found = false;
                for j in i + 1 .. placements.len(){
                    if problem.musicians[j] != tj {
                        continue;
                    }


                    scoreBefore = scoreBefore + eval_placement(problem, placements, i, false, cache) + eval_placement(problem, placements, j, false, cache);
                    let scorei = cache.placement_score[i];
                    let scorej = cache.placement_score[j];
                    swapHist.push((i, j, scorei, scorej));

                    placements.swap(i, j);
                    cache.swap_placement(i, j);
                    found = true;
                    scoreAfter =  scoreAfter +  eval_placement(problem, placements, i, false, cache) + eval_placement(problem, placements, j, false, cache);
                    break;
                }
                if found {
                    continue;
                }
                for j in i + 1 .. placements.len(){
                    if problem.musicians[j] == ti {
                        continue;
                    }
                    scoreBefore = scoreBefore + eval_placement(problem, placements, i, false, cache) + eval_placement(problem, placements, j, false, cache);
                    let scorei = cache.placement_score[i];
                    let scorej = cache.placement_score[j];
                    swapHist.push((i, j, scorei, scorej));

                    placements.swap(i, j);
                    cache.swap_placement(i, j);
                    found = true;
                    scoreAfter =  scoreAfter +  eval_placement(problem, placements, i, false, cache) + eval_placement(problem, placements, j, false, cache);
                    break;
                }
            }
            if beforeScoreStrict > eval(problem, placements, cache) {
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

fn try_swap<R:Rng>(problem:&Problem, placements:&mut Vec<Point>, rng:&mut R, cache:&mut CacheState) {
    for i in 0..placements.len()
    {
        for j in i + 1..placements.len()
        {
            if problem.musicians[i] == problem.musicians[j] {
                continue;
            }
                let score1 = eval_placement(problem, placements, i, false, cache) + eval_placement(problem, placements, j, false, cache);
                let scorei = cache.placement_score[i];
                let scorej = cache.placement_score[j];

                placements.swap(i, j);
                cache.swap_placement(i, j);
                
                let score2 = eval_placement(problem, placements, i, false, cache) + eval_placement(problem, placements, j, false, cache);
                if score1 > score2 {
                    placements.swap(i, j);
                    cache.swap_placement(i, j);
                    cache.placement_score[i] = scorei;
                    cache.placement_score[j] = scorej;
            }
        }
    } 
}

fn eval(problem:&Problem, placements:&Vec<Point>, cache:&mut CacheState) -> f64 {
    let mut result = 0.0;
    for i in 0..placements.len()
    {
        result += eval_placement(problem, placements, i, true, cache);
    }
    result
}

// 各ミュージシャンごとの観客の評価値の合算
fn eval_placement(problem:&Problem, placements:&Vec<Point>, index:usize, total:bool, cache:&mut CacheState) -> f64 {
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
    result
}

