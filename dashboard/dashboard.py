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

import drawsvg as dw
import flask
from flask import Flask, request, render_template, jsonify
from flask_cors import CORS
from werkzeug.serving import run_simple
from werkzeug.middleware.dispatcher import DispatcherMiddleware
from sqlalchemy import create_engine, VARCHAR, select
from sqlalchemy import Column, Integer, String, Float, DateTime
from sqlalchemy.orm import scoped_session, sessionmaker, declarative_base
from PIL import Image

visualizer_url = "http://35.221.99.118/repo/visualizer"
static_path = pathlib.Path(__file__).resolve().parent / 'static'
repo_path = pathlib.Path(__file__).resolve().parent.parent
problems_path = repo_path / "problem.json"
solutions_path = repo_path / "solutions"
app = Flask(__name__, static_folder=str(static_path), static_url_path='')
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0

engine = create_engine('mysql+pymysql://{user}:{password}@{host}/{db}?charset=utf8'.format(**{
    'host': os.environ.get('DB_HOST', 'localhost'),
    'db': os.environ.get('DB_NAME', 'main'),
    'user': os.environ.get('DB_USER', 'dp'),
    'password': os.environ.get('DB_PASSWORD', 'password'),
}))

CORS(
    app,
    supports_credentials=True
)

problem_file_cache = {}


@app.after_request
def add_header(response):
    if 'Expires' in response.headers:
        del response.headers['Expires']
    response.headers['Cache-Control'] = 'no-store'
    return response


def gen_thumbnail(src_path: pathlib.Path, dst_path: pathlib.Path):
    img = Image.open(src_path)
    img.resize((100, 100)).save(dst_path)


@app.route('/eval_solution', methods=["GET", "POST"])
def eval_solution():
    solution = request.args["solution"]
    fd, tmpfile = tempfile.mkstemp()
    print(tmpfile)
    with open(tmpfile, 'w+b') as fp:
        fp.write(solution.encode())
        fp.close()
    env = os.environ.copy()
    env["ISL_FILE"] = tmpfile
    env["PROBLEM_ID"] = request.args["problem_id"]
    cp = subprocess.run(["node_modules/.bin/ts-node", "index.ts"], capture_output=True, env=env, cwd="../eval-v2")
    print(cp.stdout.decode(), file=sys.stderr)
    print(cp.stderr.decode(), file=sys.stderr)
    lines = cp.stdout.decode().splitlines()
    if len(lines) == 0:
        return "failed"
    line = lines[-1]
    return line


def sort_problems(problems):
    reverse = False
    if request.args.get("desc"):
        reverse = True

    if request.args.get("sort-by"):
        key = request.args.get("sort-by")
        problems.sort(key=lambda x: x[key] if key in x else x["id"], reverse=reverse)
    return problems


def list_solutions():
    solutions = defaultdict(lambda: [])
    submission_results = list([os.path.relpath(x, solutions_path)
                               for x in glob.glob(str(solutions_path / "*.submission.result"))])
    for r in submission_results:
        with open(solutions_path / r) as f:
            js = json.load(f)
        if "Success" in js:
            solution_name = r[:-len(".submission.result")]
            problem_id = int(js["Success"]["submission"]["problem_id"])
            valid = False
            processing = False
            score = 0
            if "Success" in js["Success"]["submission"]["score"]:
                score = js["Success"]["submission"]["score"]["Success"]
                valid = True
            elif "Processing" in js["Success"]["submission"]["score"]:
                processing = True

            svg = solution_name + ".svg"
            if not os.path.exists(static_path / svg):
                print(r)
                print("A", js)
                p_js = get_problem_json(f"{problem_id}.json")
                print(solutions_path / solution_name)
                try:
                    with open(solutions_path / solution_name) as f:
                        s_js = json.load(f)
                        solution_svg(p_js, s_js).save_svg(static_path / svg)
                except Exception:
                    valid = False
                    svg = None

            solutions[problem_id].append({
                "name": solution_name,
                "svg": svg,
                "score": score,
                "valid": valid,
                "processing": processing,
            })

    for pid, s in solutions.items():
        s.sort(key=lambda x: x["score"], reverse=True)

    return solutions


def get_problem_json(problem_name):
    global problem_file_cache
    if problem_name not in problem_file_cache:
        with open(problems_path / problem_name) as f:
            js = json.load(f)
        problem_file_cache[problem_name] = js
    return problem_file_cache[problem_name]


@app.route('/solutions/<int:problem_id>')
def get_solutions(problem_id: int):
    js = get_problem_json(f"{problem_id}.json")
    # TODO: clean
    p = {}
    p["room_width"] = js["room_width"]
    p["room_height"] = js["room_height"]
    p["stage_width"] = js["stage_width"]
    p["stage_height"] = js["stage_height"]
    p["stage_bottom_left"] = js["stage_bottom_left"]
    p["musicians_size"] = len(js["musicians"])
    p["attendees_size"] = len(js["attendees"])
    p["svg"] = f"{problem_id}.svg"
    if not os.path.exists(static_path / p["svg"]):
        problem_svg(js).save_svg(static_path / p["svg"])

    solutions = list_solutions()[problem_id]
    return render_template(
        'solutions.jinja2',
        problem_id=problem_id,
        problem=p,
        solutions=solutions
    )


def problem_svg(js):
    d = dw.Drawing(js["room_width"], js["room_height"], id_prefix='id', transform='scale(1,-1)')
    d.append(dw.Rectangle(0, 0, js["room_width"], js["room_height"], fill="silver"))
    d.append(dw.Rectangle(js["stage_bottom_left"][0], js["stage_bottom_left"][1],
                          js["stage_width"], js["stage_height"], fill='#D0E0F0'))
    for a in js["attendees"]:
        d.append(dw.Circle(a["x"], a["y"], 2, fill="black"))
    return d


def solution_svg(p_js, s_js):
    d = dw.Drawing(p_js["stage_width"], p_js["stage_height"], id_prefix='id', transform='scale(1,-1)')
    d.append(dw.Rectangle(0, 0, p_js["stage_width"], p_js["stage_height"], fill='#D0E0F0'))
    for a in s_js["placements"]:
        d.append(dw.Circle(a["x"] - p_js["stage_bottom_left"][0],
                           a["y"] - p_js["stage_bottom_left"][1],
                           10, stroke="red", fill="red", fill_opacity="0.2", stroke_opacity="0.5"))
    return d


@app.route('/')
def get_index():
    problem_files = list([os.path.relpath(x, problems_path) for x in glob.glob(str(problems_path / "*.json"))])
    problem_files.sort(key=lambda x: int(x[:-5]))
    problems = [{"id": int(x[:-5]), "name": x} for x in problem_files]
    problems_dict = {int(x["id"]): x for x in problems}

    for p in problems:
        js = get_problem_json(p["name"])
        p["room_width"] = js["room_width"]
        p["room_height"] = js["room_height"]
        p["stage_width"] = js["stage_width"]
        p["stage_height"] = js["stage_height"]
        p["stage_bottom_left"] = js["stage_bottom_left"]
        p["musicians_size"] = len(js["musicians"])
        p["attendees_size"] = len(js["attendees"])
        p["svg"] = str(p["id"]) + ".svg"
        if not os.path.exists(static_path / p["svg"]):
            problem_svg(js).save_svg(static_path / p["svg"])

    solutions = list_solutions()

    sort_problems(problems)

    return render_template(
        'index.jinja2',
        is_search=request.args.get("search"),
        solutions=solutions,
        problems=problems,
        result_by_api={},
        problems_dict=problems_dict
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
