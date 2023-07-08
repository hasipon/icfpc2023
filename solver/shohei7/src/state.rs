pub struct CacheState {
    pub attendee_visiblity:Vec<i32>,
	pub placement_score:Vec<f64>,
}

impl CacheState {
    pub fn new(problem:&Problem) -> Self {
        CacheState {

        }
    }
    pub fn clear(&mut self) {
        for pointer in &mut attendee_visiblity {
            *pointer = 0;
        }
    }
}
