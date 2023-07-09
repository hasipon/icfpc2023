import csv
import glob
import math
import os
import subprocess
import pathlib
import json
from collections import defaultdict
from typing import *
from dataclasses import dataclass

import drawsvg as dw
from flask import Flask, request, render_template, jsonify
from flask_cors import CORS

visualizer_url = "http://35.221.99.118/repo/visualizer"
static_path = pathlib.Path(__file__).resolve().parent / 'static'
repo_path = pathlib.Path(__file__).resolve().parent.parent
ideal_tsv_path = repo_path / "ideal.tsv"
heatmap_path = repo_path / "heatmap"
problems_path = repo_path / "problem.json"
solutions_path = repo_path / "solutions"
app = Flask(__name__, static_folder=str(static_path), static_url_path='')
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0
CORS(app, supports_credentials=True)

problem_file_cache = {}


@app.after_request
def add_header(response):
    if 'Expires' in response.headers:
        del response.headers['Expires']
    response.headers['Cache-Control'] = 'no-store'
    return response


@dataclass
class Problem:
    id: int
    json: Any
    best_score: float = 0
    ideal_score: float = 0

    def __init__(self, problem_id: int):
        self.id = problem_id
        self.json = get_problem_json(self.id)

    @property
    def name(self) -> str:
        return f"{self.id}.json"

    @property
    def room_width(self) -> float:
        return self.json["room_width"]

    @property
    def room_height(self) -> float:
        return self.json["room_height"]

    @property
    def stage_width(self) -> float:
        return self.json["stage_width"]

    @property
    def stage_height(self) -> float:
        return self.json["stage_height"]

    @property
    def stage_bottom_left(self) -> List[float]:
        return self.json["stage_bottom_left"]

    @property
    def musicians_size(self) -> int:
        return len(self.json["musicians"])

    @property
    def attendees_size(self) -> int:
        return len(self.json["attendees"])

    @property
    def pillars_size(self) -> float:
        return len(self.json["pillars"])

    @property
    def tastes_size(self) -> int:
        return len(self.json["attendees"][0]["tastes"])

    @property
    def musicians(self) -> List:
        return self.json["musicians"]

    @property
    def attendees(self) -> List:
        return self.json["attendees"]

    @property
    def svg_path(self) -> str:
        p = prepare_problem_svg(self.id)
        return "/" + p if p else ""

    @property
    def heatmap_image_path(self) -> str:
        p = prepare_heatmap_image(self.id)
        return "/" + p if p else ""

    @property
    def diff_ideal_score(self) -> float:
        return max(0.0, self.ideal_score - self.best_score)


@dataclass
class Solution:
    name: str  # 1-sample.json
    problem_id: int
    submission_id: str
    solution_json: Any
    score: float = 0
    valid: bool = False
    processing: bool = False

    def __init__(self, result_file_path: str):
        file_name = os.path.basename(result_file_path)
        self.name = file_name[:-len(".submission.result")]
        self.problem_id = int(file_name.split("-", 1)[0])

        with open(solutions_path / self.name) as f:
            self.solution_json = json.load(f)

        with open(result_file_path) as f:
            js = json.load(f)

        if "Success" not in js:
            return

        if "Success" in js["Success"]["submission"]["score"]:
            self.score = js["Success"]["submission"]["score"]["Success"]
            self.valid = True
        elif "Processing" in js["Success"]["submission"]["score"]:
            self.processing = True

    @property
    def svg_path(self) -> str:
        return "/" + prepare_solution_svg(self.problem_id, self.name)

    @property
    def placements(self) -> List:
        return self.solution_json["placements"]


def sort_problems(problems: Dict[int, Problem], solutions: DefaultDict[int, List[Solution]]):
    reverse = False
    if request.args.get("desc"):
        reverse = True

    if request.args.get("sort-by"):
        key = request.args.get("sort-by")
        print(key)
        return dict(sorted(problems.items(), key=lambda x: getattr(x[1], key), reverse=reverse))

    return problems


def list_solutions() -> DefaultDict[int, List[Solution]]:
    solutions: DefaultDict[int, List[Solution]] = defaultdict(lambda: [])

    for r in glob.glob(str(solutions_path / "*.submission.result")):
        submission = Solution(r)
        solutions[submission.problem_id].append(Solution(r))

    for s in solutions.values():
        s.sort(key=lambda x: x.score, reverse=True)

    return solutions


def list_problems() -> Dict[int, Problem]:
    problem_ids = [int(os.path.basename(f)[:-len(".json")]) for f in glob.glob(str(problems_path / "*.json"))]
    return {pid: Problem(pid) for pid in sorted(problem_ids)}


def get_problem_json(problem_id: int):
    global problem_file_cache
    if problem_id not in problem_file_cache:
        with open(problems_path / f"{problem_id}.json") as f:
            js = json.load(f)
        problem_file_cache[problem_id] = js
    return problem_file_cache[problem_id]


def gen_problem_svg(js):
    d = dw.Drawing(js["room_width"], js["room_height"], id_prefix='id', transform='scale(1,-1)')
    d.append(dw.Rectangle(0, 0, js["room_width"], js["room_height"], fill="silver"))
    d.append(dw.Rectangle(js["stage_bottom_left"][0], js["stage_bottom_left"][1],
                          js["stage_width"], js["stage_height"], fill='#D0E0F0'))
    for p in js["pillars"]:
        d.append(dw.Circle(p["center"][0], p["center"][1], p["radius"], fill="maroon"))
    for a in js["attendees"]:
        d.append(dw.Circle(a["x"], a["y"], 2, fill="black"))
    return d


def prepare_problem_svg(problem_id):
    svg = str(problem_id) + ".svg"
    if not os.path.exists(static_path / svg):
        p_js = get_problem_json(problem_id)
        gen_problem_svg(p_js).save_svg(static_path / svg)
    return svg


def gen_solution_svg(p_js, s_js):
    d = dw.Drawing(p_js["stage_width"], p_js["stage_height"], id_prefix='id', transform='scale(1,-1)')
    d.append(dw.Rectangle(0, 0, p_js["stage_width"], p_js["stage_height"], fill='#000000'))
    d.append(dw.Rectangle(5, 5, p_js["stage_width"] - 10, p_js["stage_height"] - 10, fill='#D0E0F0'))
    for a in s_js["placements"]:
        d.append(dw.Circle(a["x"] - p_js["stage_bottom_left"][0],
                           a["y"] - p_js["stage_bottom_left"][1],
                           5, stroke="red", fill="red", fill_opacity="0.2", stroke_opacity="0.5"))
    return d


def prepare_solution_svg(problem_id, solution_name):
    svg = solution_name + ".svg"
    if not os.path.exists(static_path / svg):
        p_js = get_problem_json(problem_id)
        try:
            with open(solutions_path / solution_name) as f:
                s_js = json.load(f)
                gen_solution_svg(p_js, s_js).save_svg(static_path / svg)
        except Exception as e:
            print(e)
    return svg


def gen_heatmap_image(p_js, heatmap):
    g_stage = dw.Group(id="stage", fill='none', transform='scale(1,-1)')
    g_stage.append(dw.Rectangle(0, 0, p_js["stage_width"], p_js["stage_height"], fill='white'))
    min_score = max(x[0] for x in heatmap)
    max_score = min(x[0] for x in heatmap)

    for e in heatmap:
        score, x, y = e
        h = ((score - min_score) / (max_score - min_score)) * 240
        g_stage.append(dw.Rectangle(x - p_js["stage_bottom_left"][0],
                                    y - p_js["stage_bottom_left"][1],
                                    10, 10, fill=f'hsl({h}, 100%, 50%)'))

    wh = max(p_js["stage_width"], p_js["stage_height"]) * 0.1
    g_scale = dw.Group(id="scale", fill='none')
    for i in range(10):
        g_scale.append(dw.Rectangle(i * wh, 0, wh, wh, fill=f'hsl({(i/9) * 240}, 100%, 50%)'))
        sample_score = min_score + (i * (max_score - min_score)) / 9
        print(sample_score)
        g_scale.append(dw.Text(f"{sample_score:.2f}", wh / 5, wh / 2 + i * wh, wh / 2, center=True, fill='black'))

    wh = max(p_js["stage_width"], p_js["stage_height"])
    d = dw.Drawing(wh, wh * 1.05 + wh, id_prefix='id')
    d.append(dw.Use(g_stage, 0, p_js["stage_height"]))
    d.append(dw.Use(g_scale, 0, p_js["stage_height"] * 1.05))
    return d


def prepare_heatmap_image(problem_id):
    file_name = f"{problem_id}-heatmap.svg"
    if os.path.exists(static_path / file_name):
        return file_name
    try:
        if not os.path.exists(heatmap_path / f"{problem_id}.csv"):
            return None
        with open(heatmap_path / f"{problem_id}.csv") as f:
            heatmap = [(float(row[0]), float(row[1]), float(row[2])) for row in csv.reader(f) if len(row) == 3]
        if not heatmap:
            return None
        p_js = get_problem_json(problem_id)
        gen_heatmap_image(p_js, heatmap).save_svg(str(static_path / file_name))
    except Exception as e:
        print(e)
        return None
    return file_name


def get_ideal_scores() -> Dict[int, float]:
    ideal_scores: Dict[int, float] = dict()
    with open(ideal_tsv_path) as f:
        for row in csv.reader(f, delimiter="\t"):
            if len(row) == 2:
                ideal_scores[int(row[0])] = float(row[1])
    return ideal_scores


@app.route('/')
def get_index():
    problems = list_problems()
    solutions = list_solutions()

    for pid, prob_solution in solutions.items():
        if 0 < len(prob_solution):
            problems[pid].best_score = prob_solution[0].score

    for pid, score in get_ideal_scores().items():
        problems[pid].ideal_score = score

    problems = sort_problems(problems, solutions)

    return render_template(
        'index.jinja2',
        is_search=request.args.get("search"),
        solutions=solutions,
        problems=problems,
    )


@app.route('/solutions/<int:problem_id>')
def get_solutions(problem_id: int):
    return render_template(
        'solutions.jinja2',
        problem=Problem(problem_id),
        solutions=list_solutions()[problem_id]
    )


@app.route('/filter')
def get_filter():
    return render_template('filter.jinja2')


@app.route('/git_status')
def git_status():
    output = ""
    try:
        output += subprocess.check_output(["git", "status"], stderr=subprocess.STDOUT).decode('utf-8').strip()
    except subprocess.CalledProcessError as e:
        output += "Error:" + str(e)
    return render_template('output.jinja2', output=output)


@app.route('/git_pull')
def git_pull():
    output = ""
    try:
        output += subprocess.check_output(["git", "pull"], stderr=subprocess.STDOUT).decode(
            'utf-8').strip()
    except subprocess.CalledProcessError as e:
        output += "Error:" + str(e)
    return render_template('output.jinja2', output=output)


if __name__ == "__main__":
    app.run('0.0.0.0', port=5000, threaded=True, debug=True)
