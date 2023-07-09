
extern crate serde_json;
extern crate serde;
extern crate rand;
extern crate chrono;

mod data;
mod s_state;
mod g_state;
use data::*;
use s_state::*;
use g_state::*;
use std::{fs, env};
use rand::Rng;
use chrono::prelude::*;


fn main() -> Result<(), Box<dyn std::error::Error>> {
    let timestamp = Utc::now().timestamp();

    let args: Vec<String> = env::args().collect();
    let id = if args.len() <= 1 { "11" } else { &args[1] };
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
    let mut rng: rand::rngs::StdRng = rand::SeedableRng::seed_from_u64(3);
    
    let mut volumes = vec![1.0; problem.musicians.len()];
    let iindex:i64 = index.parse()?;
    if iindex > 55 {
        problem.extention = Option::Some(());
    }
 
    // 貪欲に配置する
    let mut placements = Vec::new();

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
    println!("{}", placements.len());
    let mut swap_state = SwapState::new(&problem);
    for i in 0..4
    {
        println!("{}: s{}:{}", index, i, s_eval(&problem, &placements, &mut volumes, &mut swap_state));
        try_swap(&problem, &mut placements, &mut rng, &mut swap_state);
        swap_state = SwapState::new(&problem);
    }

    let score = s_eval(&problem, &placements, &mut volumes, &mut swap_state);
    println!("{}:{}", index, score);

    let answer:Answer = Answer { placements, volumes };
    let answer_string = serde_json::to_string(&answer)?;
    let name = "shohei9-2";
    fs::write(
        format!("../../solutions/{}-{}.json", index, name), 
        &answer_string
    )?;
    fs::write(
        format!("../../solutions/{}-{}.myeval", index, name), 
        &score.to_string()
    )?;
    Ok(())
}
