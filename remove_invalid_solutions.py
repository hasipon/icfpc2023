import os
import glob
import json


def load_json(s):
    try:
        return json.loads(s)
    except ValueError as err:
        return False


def remove_solution(solution):
    # os.remove(solution)
    if os.path.exists(solution + ".submission"):
        os.remove(solution + ".submission")
    if os.path.exists(solution + ".submission.result"):
        os.remove(solution + ".submission.result")


def main():
    for solution in glob.glob('solutions/[1-9]*-*.json'):
        with open(solution) as f:
            s = f.read()
            js = load_json(s)
            if not js:
                print("found invalid json:", solution, "content:", s)
                remove_solution(solution)
                continue

            if "placements" not in js:
                print("found invalid json:", solution, "content:", s)
                continue

            for p in js["placements"]:
                if "x" not in p or "y" not in p:
                    print("found invalid json:", solution, "content:", s)
                    continue

                elif p["x"] is None or p["y"] is None:
                    print("found invalid json:", solution, "content:", s)
                    continue


if __name__ == "__main__":
    main()
