import glob
import os
import subprocess
import pathlib
import shutil
import math
import json
import sys
import tempfile
import urllib
from collections import defaultdict
from typing import *
from dataclasses import dataclass

import drawsvg as dw
import flask
from flask import Flask, request, render_template, jsonify
from flask_cors import CORS
from sqlalchemy import create_engine, VARCHAR, select
from sqlalchemy import Column, Integer, String, Float, DateTime
from sqlalchemy.orm import scoped_session, sessionmaker, declarative_base

visualizer_url = "http://35.221.99.118/repo/visualizer"
static_path = pathlib.Path(__file__).resolve().parent / 'static'
repo_path = pathlib.Path(__file__).resolve().parent.parent
problems_path = repo_path / "problem.json"
solutions_path = repo_path / "solutions"
app = Flask(__name__, static_folder=str(static_path), static_url_path='')
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0
CORS(app, supports_credentials=True)

engine = create_engine('mysql+pymysql://{user}:{password}@{host}/{db}?charset=utf8'.format(**{
    'host': os.environ.get('DB_HOST', 'localhost'),
    'db': os.environ.get('DB_NAME', 'main'),
    'user': os.environ.get('DB_USER', 'dp'),
    'password': os.environ.get('DB_PASSWORD', 'password'),
}))

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
    def musicians(self) -> List:
        return self.json["musicians"]

    @property
    def attendees(self) -> List:
        return self.json["attendees"]

    @property
    def svg_path(self) -> str:
        return "/" + prepare_problem_svg(self.id)


@dataclass
class Submission:
    name: str  # 1-sample.json
    problem_id: int
    submission_id: str
    solution_json: Any
    score: float = 0
    valid: bool = False
    processing: bool = False

    def __init__(self, file_path: str):
        file_name = os.path.basename(file_path)
        self.name = file_name[:-len(".submission.result")]
        self.problem_id = int(file_name.split("-", 1)[0])

        with open(solutions_path / self.name) as f:
            self.solution_json = json.load(f)

        with open(file_path) as f:
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


def sort_problems(problems):
    reverse = False
    if request.args.get("desc"):
        reverse = True

    if request.args.get("sort-by"):
        key = request.args.get("sort-by")
        problems.sort(key=lambda x: x[key] if key in x else x["id"], reverse=reverse)
    return problems


def list_solutions() -> DefaultDict[int, List[Submission]]:
    solutions: DefaultDict[int, List[Submission]] = defaultdict(lambda: [])

    for r in glob.glob(str(solutions_path / "*.submission.result")):
        submission = Submission(r)
        solutions[submission.problem_id].append(Submission(r))

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
    d.append(dw.Rectangle(0, 0, p_js["stage_width"], p_js["stage_height"], fill='#D0E0F0'))
    for a in s_js["placements"]:
        d.append(dw.Circle(a["x"] - p_js["stage_bottom_left"][0],
                           a["y"] - p_js["stage_bottom_left"][1],
                           10, stroke="red", fill="red", fill_opacity="0.2", stroke_opacity="0.5"))
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


@app.route('/')
def get_index():
    problems = list_problems()
    solutions = list_solutions()
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


@app.route('/vis/<solution>')
def get_vis(solution: str):
    problem_id, isl = engine.execute("SELECT problem_id, isl FROM solution WHERE id=%s", (solution,)).fetchone()
    return flask.redirect(visualizer_url + f"/#{problem_id};{urllib.parse.quote(isl)}")


@app.route('/eval_output/<solution>')
def eval_output(solution: str):
    (output,) = engine.execute("SELECT eval_output FROM solution WHERE id=%s", (solution,)).fetchone()
    response = flask.make_response(output, 200)
    response.mimetype = "text/plain"
    return response


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
