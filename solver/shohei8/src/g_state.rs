use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::collections::{HashMap, HashSet};

use rand::Rng;

use crate::data::*;
use crate::s_state::*;

pub struct GridState {
    pub musician_groups:HashMap<usize, Vec<usize>>,
    pub grid_sight:Vec<Vec<Sight>>, // 各座標から柱に邪魔されない客
    pub x:f64,
    pub y:f64,
    pub w:i32,
    pub h:i32,
    pub scale_x:f64,
    pub scale_y:f64,
}

pub struct GridCache {
    pub filled  :HashMap<i32, usize>,
    pub crossing:HashMap<i32, Vec<Pair>>,
}

pub struct Pair {
    m:usize,
    a:usize,
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct P {
    x:i32,
    y:i32,
}

impl GridState {
    pub fn new(problem:&Problem) -> GridState {
        let x = problem.stage_bottom_left.0 + 5.0;
        let y = problem.stage_bottom_left.1 + 5.0;
        let width = problem.stage_width;
        let height = problem.stage_height;
        let w = (width  / 10.0).floor() as i32;
        let h = (height / 10.0).floor() as i32;
        let scale_x = (width  - 10.0) / (w - 1) as f64; 
        let scale_y = (height - 10.0) / (h - 1) as f64; 
        let mut musician_groups:HashMap<usize, Vec<usize>> = HashMap::new();
        
        for i in 0..problem.musicians.len() 
        {
            let group = musician_groups.get(&problem.musicians[i]);
            if group.is_none() {
                musician_groups.insert(i, Vec::new());
            }
            musician_groups.get_mut(&problem.musicians[i]).unwrap().push(i);
        }
        
        let grid_sight = Vec::new();
        let mut result = GridState { 
            musician_groups,
            grid_sight,
            w,
            h,
            scale_x,
            scale_y,
            x,
            y
        };
        for iy in 0..h 
        {
            for ix in 0..w 
            {
                let p = P{x:ix, y:iy};
                result.grid_sight.push(resolve_sight(problem, &SightMode::Init, p.to_point(&result)));
            }
        }
        result        
    } 

    pub fn init_grid(&self, problem:&Problem) -> Vec<Point> {
        let mut result: Vec<Point> = Vec::new();
        let mut cache = GridCache::new(self, problem);
        for i in 0..problem.musicians.len()
        {
            let best = self.find_best(problem, &result, i, &mut cache);
            cache.add(&self, &problem, best, i);
            result.push(best.to_point(&self));
        }
        result
    }
    
    pub fn find_best(&self, problem:&Problem, placements:&Vec<Point>, index:usize, cache:&mut GridCache) -> P {
        let mut best = P {x:0, y:0};
        let mut best_score = -std::f64::INFINITY;
        println!("{}/{}", index, problem.musicians.len());
        for x in 0..self.w {
            for y in 0..self.h {
                let p = P {x, y};
                let center = p.to_point(self);
                if cache.filled.contains_key(&p.to_index(self)) {
                    continue;
                }
                let mut score = 0.0;
                let sights = resolve_sight(problem, &SightMode::Placements(placements), center);

                for sight in sights {
                    let a = &problem.attendees[sight.attendee];
                    score += 1000000.0 * a.tastes[problem.musicians[index]] / sight.d2;
                }
                if problem.extention.is_some() {
                    let mut q = 1.0;

                    for i in self.musician_groups.get(&problem.musicians[index]).unwrap()
                    {
                        if *i == index { continue; }
                        let p = placements[*i];
                        let dx = center.x - p.x;
                        let dy = center.y - p.y;
                        let d = (dx * dx + dy * dy).sqrt();
                        let value = 1.0 / d.max(10.0);
                        q += value;
                    }
                    score *= q;
                }
                if score > best_score {
                    best_score = score;
                    best = p;
                }
            }
        }
        best
    }

    pub fn init_random_grid<R:Rng>(&mut self, problem:&Problem, rng:&mut R) -> Vec<Point> {
        let mut set = HashSet::new();
        let mut result: Vec<Point> = Vec::new();
        for _ in &problem.musicians
        {
            loop {
                let x = rng.gen_range(0..self.w);
                let y = rng.gen_range(0..self.h);
                let key = x + y * self.w;
                if !set.contains(&key)
                {
                    set.insert(key);
                    result.push(Point { 
                        x: self.x + (x as f64 * self.scale_x) + 0.5, 
                        y: self.y + (y as f64 * self.scale_y) + 0.5,
                    });
                    break;
                }
            }
        }
        result
    }

    fn to_point(&self, placement:&Point) -> P {
        P {
            x:((placement.x - self.x) / self.scale_x).floor() as i32,
            y:((placement.y - self.y) / self.scale_y).floor() as i32,
        }
    }
}

impl GridCache {
    pub fn new(state:&GridState, problem:&Problem) -> GridCache {
        GridCache {
            filled: HashMap::new(), 
            crossing: HashMap::new()
        }
    }

    pub fn add(&mut self, state:&GridState, problem:&Problem, point:P, index:usize) {
        self.filled.insert(point.to_index(state), index);

        // TODO: crossing
    }
}

impl P {
    fn to_point(&self, state:&GridState) -> Point {
        Point {
            x: state.x + self.x as f64 * state.scale_x,
            y: state.y + self.y as f64 * state.scale_y,
        }
    }
    fn to_index(&self, state:&GridState) -> i32 {
        return self.x + self.y * state.w
    }
}

enum SightMode<'a> {
    Init,
    Placements(&'a Vec<Point>),
}

fn resolve_sight(problem:&Problem, mode:&SightMode, center:Point) -> Vec<Sight> {
    let mut nearest_d = std::f64::INFINITY;
    let mut nearest_dir = 0.0;
    let mut nodes = Vec::new();
    
    match mode {
        &SightMode::Init => {},
        &SightMode::Placements(placements) => {
            
            // ミュージシャンのふちを追加
            for (i, p) in placements.iter().enumerate()
            {
                if *p == center { continue; }
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
        },
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
    sights
}