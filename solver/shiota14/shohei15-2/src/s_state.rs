use std::{vec, cmp::Ordering, collections::HashMap};

use crate::data::*;
use std::collections::BinaryHeap;
use rand::Rng;

// スワップ処理用の実装
pub struct SwapState {
    pub musician_groups:HashMap<usize, Vec<usize>>,
	pub placement_score:Vec<f64>,
    pub placement_sights:Vec<Option<Vec<Sight>>>,
}

impl SwapState {
    pub fn new(problem:&Problem) -> Self {
        let mut placement_sights = Vec::new();
        let mut musician_groups:HashMap<usize, Vec<usize>> = HashMap::new();
        for i in 0..problem.musicians.len() 
        {
            placement_sights.push(Option::None);
            let group = musician_groups.get(&problem.musicians[i]);
            if group.is_none() {
                musician_groups.insert(i, Vec::new());
            }
            musician_groups.get_mut(&problem.musicians[i]).unwrap().push(i);
        }
        SwapState {
            placement_score: vec![0.0; problem.musicians.len()],
            placement_sights,
            musician_groups
        }
    }

    pub fn swap_placement(&mut self, a:usize, b:usize) {
        self.placement_score .swap(a, b);
        self.placement_sights.swap(a, b);
    }
}


#[derive(Clone, Copy, PartialEq, Debug)]
pub struct EvalNode {
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

#[derive(Clone, Copy, PartialOrd, PartialEq, Debug)]
pub enum EvalNodeKind{
    Start(f64),
    Center(usize),
}
pub struct Sight {
    pub attendee:usize,
    pub d2:f64,
}

pub fn s_eval(problem:&Problem, placements:&Vec<Point>, volumes:&mut Vec<f64>) -> f64 {
    let mut cache = SwapState::new(&problem);
    let mut result = 0.0;
    for i in 0..placements.len()
    {
        let (score, volume) = s_eval_placement(problem, placements, i, true, &mut cache);
        result += score;
        volumes[i] = volume;
    }
    result
}

// 各ミュージシャンごとの観客の評価値の合算
pub fn s_eval_placement(problem:&Problem, placements:&Vec<Point>, index:usize, total:bool, cache:&mut SwapState) -> (f64, f64) {
    let mut nearest_d = std::f64::INFINITY;
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
        result += (1000000.0 * a.tastes[problem.musicians[index]] / sight.d2).ceil();
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
    if result < 0.0 {(0.0, 0.0)} else {(result * 10.0, 10.0)}
}


pub fn try_swap<R:Rng>(problem:&Problem, placements:&mut Vec<Point>, rng:&mut R, cache:&mut SwapState) {
    
    for i in 0..placements.len()
    {
        for j in i + 1..placements.len()
        {
            if true { 
                let score1 = s_eval_placement(problem, placements, i, false, cache).0 + s_eval_placement(problem, placements, j, false, cache).0;
                let scorei = cache.placement_score[i];
                let scorej = cache.placement_score[j];

                placements.swap(i, j);
                cache.swap_placement(i, j);
                
                let score2 = s_eval_placement(problem, placements, i, false, cache).0 + s_eval_placement(problem, placements, j, false, cache).0;
                if score1 > score2 {
                    placements.swap(i, j);
                    cache.swap_placement(i, j);
                    cache.placement_score[i] = scorei;
                    cache.placement_score[j] = scorej;
                }
            }
        }
    } 
}
