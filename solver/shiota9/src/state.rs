use std::{vec, cmp::Ordering, collections::HashMap};

use crate::data::Problem;

pub struct CacheState {
    pub musician_groups:HashMap<usize, Vec<usize>>,
	pub placement_score:Vec<f64>,
    pub placement_sights:Vec<Option<Vec<Sight>>>,
}

impl CacheState {
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
        // println!("{:?}", musician_groups);
        CacheState {
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