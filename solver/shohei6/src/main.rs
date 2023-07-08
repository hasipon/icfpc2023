
extern crate serde_json;
extern crate serde;
extern crate rand;
extern crate chrono;

mod data;
use data::*;
use std::{fs, cmp::Ordering, env};
use rand::Rng;
use chrono::prelude::*;
use std::collections::BinaryHeap;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let timestamp = Utc::now().timestamp();

    let args: Vec<String> = env::args().collect();
    let id = if args.len() <= 1 { "81" } else { &args[1] };
    solve(id, timestamp)?;

    Ok(())
}

fn solve(index:&str, timestamp:i64) -> Result<(), Box<dyn std::error::Error>> {
    let data:String = fs::read_to_string(format!("../../problem.json/{}.json", index))?; 
    let mut problem:Problem = serde_json::from_str(&data)?;
    let mut placements = Vec::new();
    let x = problem.stage_bottom_left.0 + 10.0;  
    let y = problem.stage_bottom_left.1 + 10.0;
    let w = problem.stage_width - 20.0;
    let h = problem.stage_height - 20.0;
    let mut max_taste = 0.0;
    let mut rng: rand::rngs::StdRng = rand::SeedableRng::seed_from_u64(3);
    let iindex:i64 = index.parse()?;
    if iindex > 55 {
        problem.extention = Option::Some(());
    }

    // 外側に配置する
    {
        let mut dx = 0.0;
        let mut dy: f64 = 1.0;
        let mut left = x;
        let mut right = x + w;
        let mut top = y;
        let mut bottom = y + h;
        let mut cx = right;
        let mut cy = y;

        let mut i = 0;
        while i < problem.musicians.len()
        {
            let speed = 10.0;
            placements.push(Point{x:cx, y:cy});
            i += 1;
            
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
    try_swap(&problem, &mut placements, &mut rng);
    println!("{}:{}", index, eval(&problem, &placements));

    let answer:Answer = Answer { placements };
    let answer_string = serde_json::to_string(&answer)?;
    fs::write(
        format!("../../solutions/{}-shohei6-3.json", index), 
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
        placements[i] = p1;
    }
    !hit
}

fn yama<R:Rng>(
    x:f64,
    y:f64,
    w:f64,
    h:f64,
    problem:&Problem, placements:&mut Vec<Point>, 
    power:f64,
    rng:&mut R) {
    for i in 0..placements.len()
    {
        yama_placement(x, y, w, h, problem, placements, i, power, rng);
    } 
}

fn try_swap<R:Rng>(problem:&Problem, placements:&mut Vec<Point>, rng:&mut R) {
    let rate = (500.0 / placements.len() as f64).min(1.0);
    for i in 0..placements.len()
    {
        println!("{}/{}", i, placements.len());
        for j in i + 1..placements.len()
        {
            let rand = rng.gen_range(0.0..1.0); 
            if true { 
                let score1 = eval_placement(problem, placements, i, false) + eval_placement(problem, placements, j, false);
                placements.swap(i, j);

                let score2 = eval_placement(problem, placements, i, false) + eval_placement(problem, placements, j, false);
                if score1 > score2 {
                    placements.swap(i, j);
                }
            }
        }
    } 
}
fn yama_placement<R:Rng>(
    x:f64,
    y:f64,
    w:f64,
    h:f64,
    problem:&Problem, 
    placements:&mut Vec<Point>, 
    index:usize, 
    power:f64,
    rng:&mut R) {
    let center = placements[index];
    let mut result = center;
    let mut max = eval_placement(problem, placements, index, false);
    for offset in [
        (0.0, 0.9), 
        (0.0, -0.9), 
        (0.9, 0.0), 
        (-0.9, 0.0),
        (1.0 * rng.gen_range(0.0..1.0), 1.0 * rng.gen_range(0.0..1.0)), 
        (1.0 * rng.gen_range(0.0..1.0), -1.0 * rng.gen_range(0.0..1.0)), 
        (-1.0 * rng.gen_range(0.0..1.0), 1.0 * rng.gen_range(0.0..1.0)), 
        (-1.0 * rng.gen_range(0.0..1.0), -1.0 * rng.gen_range(0.0..1.0))
    ]
    {
        let mut px = center.x + offset.0 * power;
        let mut py = center.y + offset.1 * power;
        placements[index] = Point{
            x: px,
            y: py,
        };
        let score = eval_placement(problem, placements, index, false);
        if score > max {
            max = score;
            result = placements[index];
        }
    }

    placements[index] = result;
}

fn eval(problem:&Problem, placements:&Vec<Point>) -> f64 {
    let mut result = 0.0;
    for i in 0..placements.len()
    {
        result += eval_placement(problem, placements, i, true);
    }
    result
}

// 各ミュージシャンごとの観客の評価値の合算
fn eval_placement(problem:&Problem, placements:&Vec<Point>, index:usize, total:bool) -> f64 {
    let mut nearest_d = 10000000000000000000000.0;
    let mut nearest_dir = 0.0;
    let center = placements[index];
    let mut nodes = Vec::new();
    let mut q = 1.0;

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
        if problem.musicians[i] == problem.musicians[index] {
            q += 1.0 / d.max(10.0);
        }
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
    let mut result = 0.0;
    for node in nodes
    {
        let next_dir = node.dir;
        while !heap.peek().is_none() && heap.peek().unwrap().end_dir() < next_dir
        {
            heap.pop();
        }
        match &node.kind
        {
            &EvalNodeKind::Start(_) => { heap.push(node); }
            &EvalNodeKind::Center(aindex) => {
                let a = &problem.attendees[aindex];
                if heap.peek().is_none() || node.d < heap.peek().unwrap().d {
                    result += 1000000.0 * a.tastes[problem.musicians[index]] / node.d / node.d;
                }
            },
        }
    }
    if problem.extention.is_some() {
        result *= q;
    }
    result
}

#[derive(Clone, PartialEq, Debug)]
struct EvalNode {
    pub d:f64,
    pub dir:f64,
    pub kind:EvalNodeKind
}
impl PartialOrd for EvalNode {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        other.d.partial_cmp(&self.d)
    }
} 

impl Ord for EvalNode {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(&other).unwrap()
    }
}
impl Eq for EvalNode { }

impl EvalNode {
    pub fn end_dir(&self) ->f64
    {
        match &self.kind
        {
            &EvalNodeKind::Start(dir) => { self.dir + dir }
            &EvalNodeKind::Center(_) => { self.dir }
        }
    }
}

#[derive(Clone, PartialOrd, PartialEq, Debug)]
enum EvalNodeKind{
    Start(f64),
    Center(usize),
}
