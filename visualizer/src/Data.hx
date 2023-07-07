package ;

/**
 * ...
 * @author shohei909
 */
typedef Problem = 
{
	room_width: Float,
	room_height: Float,
	stage_width: Float,
	stage_height: Float,
	stage_bottom_left: Array<Float>,
	musicians: Array<Int>,
	attendees:Array<Attendee>,
}

typedef Attendee = {
	x:Float,
	y:Float,
	tastes:Array<Float>
}