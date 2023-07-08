package;
import haxe.DynamicAccess;
import haxe.Http;
import haxe.Json;
import js.Browser;
import js.html.CanvasElement;
import js.html.Element;
import js.html.InputElement;
import js.html.TextAreaElement;
import pixi.core.Application;
import pixi.core.graphics.Graphics;
import data.Data;

class Main
{
	static var mainCanvas:CanvasElement;
	static var mainPixi:Application;
	static var currentProblem:Int;
	static var nextProblem :Int;
	static var problemInput:InputElement;
	static var input       :TextAreaElement;
	static var problem     :Problem;
	static var problemGraphics:Graphics;
	static var answerGraphics :Graphics;
	static var scale          :Float;
	static var answer:Answer;
	
	static function main() 
	{
		currentProblem = -1;
		nextProblem = -1;
		mainCanvas = cast Browser.document.getElementById("pixi");
		mainPixi = new Application({
			view  :mainCanvas,
			transparent: true,
			width: mainCanvas.width,
			height: mainCanvas.height,
			autoResize: true,
		});
		mainPixi.stage.addChild(problemGraphics = new Graphics());
		mainPixi.stage.addChild(answerGraphics  = new Graphics());
		
		problemInput = cast Browser.document.getElementById("problem");
		problemInput.onchange = onProblemChanged;
		problemInput.oninput = onProblemChanged;
		
		input = cast Browser.document.getElementById("input");
		input.onchange = onInputChanged;
		input.oninput = onInputChanged;
		
		Browser.window.addEventListener('hashchange', readHash, false);
		readHash();
	}
	
	static function readHash():Void 
	{
		var hash = Browser.location.hash;
		var index = hash.indexOf(";");
		if (index > 0)
		{
			var inputText = hash.substr(index + 1);
			input.value = StringTools.urlDecode(inputText);
		}
		else
		{
			index = hash.length;
		}
		var problem = Std.parseInt(hash.substring(1, index));
		if (problem != null && problem > 0)
		{
			problemInput.value = "" + problem;
		}
		
		onProblemChanged();
		setHash();
	}
	
	static function onProblemChanged():Void
	{
		var nextProblem = Std.parseInt(problemInput.value);
		if (Main.nextProblem == nextProblem) { return; }
		Main.nextProblem = nextProblem;
		
		var url = "../problem.json/" + nextProblem + ".json";
		var http = new Http(url);
		http.onData = onProblemGet.bind(nextProblem);
		http.request();
	}
	
	static function onProblemGet(problemId:Int, json:String):Void
	{
		if (problemId == currentProblem) { return; }
		if (problemId != nextProblem) { return; }
		currentProblem = problemId;
		
		problem = Json.parse(json);
		var w = 1.0;
		var h = 1.0;
		for (a in problem.attendees)
		{
			if (w < a.x) { w = a.x; }
			if (h < a.y) { h = a.y; }
		}
		var scaleX = 720 / w;
		var scaleY = 480 / h;
		scale = Math.min(scaleX, scaleY);
		problemGraphics.clear();
		problemGraphics.beginFill(0x000000, 1);
		for (a in problem.attendees)
		{
			problemGraphics.drawCircle(
				a.x * scale, 
				480 - a.y * scale, 
				1
			);
		}
		problemGraphics.endFill();
		problemGraphics.beginFill(0x0000FF, 0.1);
		problemGraphics.lineStyle(1, 0x0000FF, 1);
		problemGraphics.drawRect(
			problem.stage_bottom_left[0] * scale,
			480 - (problem.stage_height + problem.stage_bottom_left[1]) * scale,
			problem.stage_width * scale,
			problem.stage_height * scale
		);
		
		onInputChanged();
	}
	
	static function onInputChanged():Void
	{
		answer = Json.parse(input.value);
		answerGraphics.clear();
		answerGraphics.lineStyle(1, 0xFF0000, 0.5);
		answerGraphics.beginFill(0xFF0000, 0.2);
		
		for (p in answer.placements)
		{
			answerGraphics.drawCircle(
				p.x * scale, 
				480 - p.y * scale, 
				5 * scale
			);
		}
		setHash();
	}
	
	private static function setHash():Void 
	{
		Browser.location.hash = "#" + Std.parseInt(problemInput.value) + ";" + StringTools.urlEncode(input.value);
	}
}
