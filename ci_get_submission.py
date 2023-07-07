import os
import glob
import json
import urllib.request
import urllib.parse

AWESOME_UA='Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/114.0.0.0 Safari/537.36'


def get_submission(submission_id: int):
    assert os.getenv("ICFPC_TOKEN")
    token = os.getenv("ICFPC_TOKEN")
    request = urllib.request.Request(
        "https://api.icfpcontest.com/submission?submission_id={}".format(submission_id),
        headers={'Authorization': 'Bearer {}'.format(token), 'User-Agent': AWESOME_UA, 'Content-Type': 'application/json'})
    response = urllib.request.urlopen(request)
    print(response.getcode())
    response_body = response.read()
    print(response_body)
    return response_body.decode('utf-8')


def main():
    submissions = set(sorted(glob.glob('solutions/[1-9]*-*.json.submission')))
    results = set(sorted(glob.glob('solutions/[1-9]*-*.json.submission.result')))
    submission_not_get_result = [x for x in submissions if x + ".result" not in results]

    print(submission_not_get_result)

    for submission in submission_not_get_result:
        with open(submission) as f:
            submission_id = f.read().strip()
        response = get_submission(submission_id)
        with open('{}.result'.format(submission), mode="w") as result_file:
            result_file.write(response.strip("\""))
            result_file.write("\n")
        exit(0)


if __name__ == "__main__":
    main()
