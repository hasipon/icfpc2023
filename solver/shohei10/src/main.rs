
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
    let id = if args.len() <= 1 { "81" } else { &args[1] };
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

    // ランダムに配置する
    let mut grid_state = GridState::new(&problem);
    let mut placements = grid_state.init_grid(&problem);

    println!("{}", placements.len());
    let mut swap_state = SwapState::new(&problem);
    for i in 0..4
    {
        println!("{}: s{}:{}", index, i, s_eval(&problem, &placements, &mut volumes, &mut swap_state));
        try_swap(&problem, &mut placements, &mut rng, &mut swap_state);
        swap_state = SwapState::new(&problem);
    
        println!("{}: g{}:{}", index, i, s_eval(&problem, &placements, &mut volumes, &mut swap_state));
        //try_grid_move(&problem, &mut placements, &mut rng);
    }

    let score = s_eval(&problem, &placements, &mut volumes, &mut swap_state);
    println!("{}:{}", index, score);

    let answer:Answer = Answer { placements, volumes };
    let answer_string = serde_json::to_string(&answer)?;
    let name = "shohei10-8";
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
