use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::collections::{HashMap, HashSet};

use rand::Rng;
use rand::seq::SliceRandom;

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
    pub crossing:Vec<HashSet<Pair>>, // 各座標をまたぐ線
    pub sights  :Vec<HashMap<usize, f64>>, // 各ミュージシャンから見える客
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
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
        let width = problem.stage_width - 10.0;
        let height = problem.stage_height - 10.0;
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
    
    pub fn init_grid<R:Rng>(&self, problem:&Problem, rng:&mut R) -> Vec<Point> {
        let mut result: Vec<Point> = vec![Point{x:std::f64::NAN, y:std::f64::NAN}; problem.musicians.len()];
        let mut cache = GridCache::new(self, problem);
        let mut indexes:Vec<usize> = (0..problem.musicians.len()).collect();
        indexes.shuffle(rng);
        for i in indexes
        {
            let best = self.find_best(problem, &result, i, &mut cache);
            cache.add(&self, &problem, &result, best.0, best.1,i);
            result[i] = best.0.to_point(&self);
        }
        result
    }
    pub fn find_best(&self, problem:&Problem, placements:&Vec<Point>, index:usize, cache:&mut GridCache) -> (P, Vec<Sight>) {
        let mut best = P {x:0, y:0};
        let mut best_score = -std::f64::INFINITY;
        let mut best_sights = Vec::new();
        
        for x in 0..self.w {
            for y in 0..self.h {
                let p = P {x, y};
                let center = p.to_point(self);
                if cache.filled.contains_key(&p.to_index(self)) {
                    continue;
                }
                let mut score = 0.0;
                let sights = resolve_sight(problem, &SightMode::Placements(placements), center);

                for sight in &sights {
                    let a = &problem.attendees[sight.attendee];
                    score += 1000000.0 * a.tastes[problem.musicians[index]] / sight.d2;
                }
                if score < 0.0 { score *= 0.0000001; }
                if problem.extention.is_some() {
                    let mut q = 1.0;

                    for i in self.musician_groups.get(&problem.musicians[index]).unwrap()
                    {
                        if *i == index { continue; }
                        if *i >= placements.len() { continue; }
                        let p = placements[*i];
                        if p.y.is_nan() { continue; }
                        let dx = center.x - p.x;
                        let dy = center.y - p.y;
                        let d = (dx * dx + dy * dy).sqrt();
                        let value = 1.0 / d.max(10.0);
                        q += value;
                    }
                    score *= q;
                }
                for crossing in &cache.crossing[p.to_index(self) as usize] {
                    if let Some(d2) =  cache.sights[crossing.m].get(&crossing.a)
                    {
                        let a = &problem.attendees[crossing.a];
                        score -= 10000.0 * a.tastes[problem.musicians[crossing.m]] / d2;
                    }
                }

                if score > best_score {
                    best_score = score;
                    best = p;
                    best_sights = sights;
                }
            }
        }
        (best, best_sights)
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
                    let p = Point { 
                        x: self.x + (x as f64 * self.scale_x) + 5.0, 
                        y: self.y + (y as f64 * self.scale_y) + 5.0,
                    };
                    result.push(p);
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
        let mut sights = Vec::new();
        for _ in &problem.musicians
        {
            sights.push(HashMap::new());
        }

        let mut crossing = Vec::new();
        for _ in 0..state.w
        {
            for _ in 0..state.h
            {
                crossing.push(HashSet::new());
            }
        }
        GridCache {
            filled  : HashMap::new(), 
            crossing,
            sights,
        }
    }

    pub fn add(&mut self, state:&GridState, problem:&Problem, placements:&Vec<Point>, point:P, sights:Vec<Sight>, index:usize) {
        // 交差済みの線を削除
        let erase_targets = std::mem::replace(&mut self.crossing[point.to_index(state) as usize], HashSet::new());
        for crossing in &erase_targets
        {
            self.update_sight(state, problem, state.to_point(&placements[crossing.m]), crossing.m, crossing.a, true);
        }
        self.filled.insert(point.to_index(state), index);
        // 削除した線を再描画
        for crossing in &erase_targets
        {
            self.update_sight(state, problem, state.to_point(&placements[crossing.m]), crossing.m, crossing.a, false);
        }
        
        for sight in sights {
            self.sights[index].insert(sight.attendee, sight.d2);
            self.update_sight(state, problem, point, index, sight.attendee, false);
        }
    }
    pub fn update_sight(&mut self, state:&GridState, problem:&Problem, point:P, index:usize, attendee:usize, is_erase:bool) {
        
        // ブレゼンハムのアルゴリズム
        let center = point.to_point(state);
        let mut x0 = (center.x - state.x) / state.scale_x;
        let mut y0 = (center.y - state.y) / state.scale_y;
        let mut x1 = (problem.attendees[attendee].x - state.x) / state.scale_x;
        let mut y1 = (problem.attendees[attendee].y - state.y) / state.scale_y;

        let steep  = (y1 - y0).abs() > (x1 - x0).abs();
        let mut w = state.w;
        let mut h = state.h;
        let mut ix = point.x;
        let mut iy = point.y;

        if steep {
            std::mem::swap(&mut x0, &mut y0);
            std::mem::swap(&mut x1, &mut y1);
            std::mem::swap(&mut ix, &mut iy);
            std::mem::swap(&mut w, &mut h);
        }
        
		let dx = (x1 - x0).abs();
        let dy = (y1 - y0).abs();
        let mut error = dx / 2.0;
        let inc = if x0 < x1 { 1 } else { -1 };
        let ystep = if y0 < y1 { 1 } else { -1 }; 
        error -= dy / 2.0;

        loop {
            if error < 0.0 {
                iy += ystep;
                if iy < 0 || iy >= h {
                    self.sights[index].insert(attendee, problem.attendees[attendee].distance2(center));
                    break;
                }
                let filled = if steep { self.update_crossing(state, problem, index, attendee, center, iy, ix, is_erase) } else { self.update_crossing(state, problem, index, attendee, center, ix, iy, is_erase) };
                if filled { 
                    self.sights[index].remove(&attendee);
                    break; 
                }
                error += dx;
            }
            
            ix += inc;
            if ix < 0 || ix >= w {
                self.sights[index].insert(attendee, problem.attendees[attendee].distance2(center));
                break;
            }

            let filled = if steep { self.update_crossing(state, problem, index, attendee, center, iy, ix, is_erase) } else { self.update_crossing(state, problem, index, attendee, center, ix, iy, is_erase) };
            if filled { 
                self.sights[index].remove(&attendee);
                break; 
            }
            error -= dy;
        }
    }
    pub fn update_crossing(&mut self, state:&GridState, problem:&Problem, index:usize, attendee:usize, start:Point, x:i32, y:i32, is_erase:bool) -> bool {
        let p = P {x, y};
        let current = p.to_point(state);
        let a = &problem.attendees[attendee];
        let dx0 = a.x - start.x;
        let dy0 = a.y - start.y;
        let dx1 = current.x - start.x;
        let dy1 = current.y - start.y;
        let d = (dx0 * dx0 + dy0 * dy0).sqrt();
        let s = (dx0 * dy1 - dx1 * dy0).abs();
        let height = s / d;
        
        let hit:bool = height < 5.0;
        if is_erase {
            self.crossing[p.to_index(state) as usize].remove(&Pair{ a: attendee, m: index });
        }
        else 
        {
            if hit 
            {
                self.crossing[p.to_index(state) as usize].insert(Pair{ a: attendee, m: index });
            }    
        }
        let result = hit && self.filled.contains_key(&p.to_index(state));
        result
    }
}

impl P {
    fn to_point(&self, state:&GridState) -> Point {
        Point {
            x: state.x + self.x as f64 * state.scale_x + 5.0,
            y: state.y + self.y as f64 * state.scale_y + 5.0,
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
                if p.x.is_nan() { continue; }
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