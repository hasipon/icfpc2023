import os
import glob
import json
import urllib.request
import urllib.parse


def submit(problem_id: int, contents: str):
    assert os.getenv("ICFPC_TOKEN")
    token = os.getenv("ICFPC_TOKEN")
    data = {
        "problem_id": problem_id,
        "contents": contents,
    }
    request = urllib.request.Request(
        "https://api.icfpcontest.com/submission",
        json.dumps(data).encode(),
        headers={'Authorization': 'Bearer {}'.format(token)})
    response = urllib.request.urlopen(request)
    print(response.getcode())
    response_body = response.read()
    print(response_body)
    return response_body


def main():
    solutions = set(sorted(glob.glob('solutions/[1-9]*-*.json')))
    submissions = set(sorted(glob.glob('solutions/[1-9]*-*.json.submission')))
    solutions_not_submitted = [x for x in solutions if x + ".submission" not in submissions]

    print(solutions_not_submitted)

    for solution in solutions_not_submitted:
        problem_id = os.path.basename(solution).split("-", 1)[0]
        with open(solution) as f:
            contents = f.read()
            response = submit(problem_id, contents)
            print(type(response))
            #TODO save response to .submission


if __name__ == "__main__":
    main()