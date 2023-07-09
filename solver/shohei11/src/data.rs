use serde::{Serialize, Deserialize};

#[derive(Serialize, Deserialize, Debug)]
pub struct Problem {
	pub room_width: f64,
	pub room_height: f64,
	pub stage_width: f64,
	pub stage_height: f64,
	pub stage_bottom_left:(f64, f64),
	pub musicians:Vec<usize>,
	pub attendees:Vec<Attendee>,
	pub pillars:Vec<Pillar>,
	pub extention:Option<()>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Pillar {
	pub center: (f64, f64),
	pub radius: f64,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Attendee {
    pub x:f64,
    pub y:f64,
    pub tastes:Vec<f64>
}

fn distance2() {
	let 
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Answer {
    pub placements:Vec<Point>,
	pub volumes:Vec<f64>,
}

#[derive(Serialize, Deserialize, Debug, Copy, Clone, PartialEq)]
pub struct Point {
    pub x:f64,
    pub y:f64,
}
