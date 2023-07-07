import os
import glob
import json
import urllib.request
import urllib.parse

AWESOME_UA='Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/114.0.0.0 Safari/537.36'


def submit(problem_id: int, contents: str):
    assert os.getenv("ICFPC_TOKEN")
    token = os.getenv("ICFPC_TOKEN")
    data = {
        "problem_id": int(problem_id),
        "contents": contents,
    }
    request = urllib.request.Request(
        "https://api.icfpcontest.com/submission",
        json.dumps(data).encode(),
        headers={'Authorization': 'Bearer {}'.format(token), 'User-Agent': AWESOME_UA, 'Content-Type': 'application/json'})
    response = urllib.request.urlopen(request)
    print(response.getcode())
    response_body = response.read()
    print(response_body)
    return response_body.decode('utf-8')


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
            with open('{}.submission'.format(solution), mode="w") as subfile:
                subfile.write(response.strip("\""))
                subfile.write("\n")


if __name__ == "__main__":
    main()
