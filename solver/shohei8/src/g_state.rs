use std::collections::{HashMap, HashSet};

use rand::Rng;

use crate::data::*;



pub struct GridState {
    pub musician_groups:HashMap<usize, Vec<usize>>,
    pub x:f64,
    pub y:f64,
    pub w:i32,
    pub h:i32,
    pub scale_x:f64,
    pub scale_y:f64,
}

impl GridState {
    pub fn new(problem:&Problem) -> GridState {
        let x = problem.stage_bottom_left.0 + 5.0;
        let y = problem.stage_bottom_left.1 + 5.0;
        let width = problem.stage_width - 10.0;
        let height = problem.stage_height - 10.0;
        let w = (width  / 10.0).floor() as i32;
        let h = (height / 10.0).floor() as i32;
        let scale_x = (width  - 10.0) / (w - 1) as f64 / 10.0; 
        let scale_y = (height - 10.0) / (h - 1) as f64 / 10.0; 
        let mut musician_groups:HashMap<usize, Vec<usize>> = HashMap::new();
        
        for i in 0..problem.musicians.len() 
        {
            let group = musician_groups.get(&problem.musicians[i]);
            if group.is_none() {
                musician_groups.insert(i, Vec::new());
            }
            musician_groups.get_mut(&problem.musicians[i]).unwrap().push(i);
        }
        
        GridState { 
            musician_groups,
            w,
            h,
            scale_x,
            scale_y,
            x,
            y
        }
    } 

    pub fn init_grid<R:Rng>(&mut self, problem:&Problem, rng:&mut R) -> Vec<Point> {
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
                        x: self.x + x as f64 * self.scale_x, 
                        y: self.y + y as f64 * self.scale_y,
                    });
                    break;
                }
            } 
        }
        result
    }
    
    pub fn try_yama_grid_move<R:Rng>(&mut self, problem:&Problem, placement:&Vec<Point>, rng:&mut R) {
    }
}


