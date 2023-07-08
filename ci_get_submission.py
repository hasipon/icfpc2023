import os
import glob
import json
import urllib.request
import urllib.parse
import time

AWESOME_UA='Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/114.0.0.0 Safari/537.36'


def get_submission(submission_id):
    print("get_submission", submission_id)
    assert os.getenv("ICFPC_TOKEN")
    token = os.getenv("ICFPC_TOKEN")
    request = urllib.request.Request(
        "https://api.icfpcontest.com/submission?submission_id={}".format(submission_id),
        headers={'Authorization': 'Bearer {}'.format(token), 'User-Agent': AWESOME_UA, 'Content-Type': 'application/json'})
    response = urllib.request.urlopen(request)
    print(response.getcode())
    response_body = response.read().decode('utf-8')
    print(response_body)
    return response_body


def main():
    submissions = set(sorted(glob.glob('solutions/[1-9]*-*.json.submission')))
    results = set(sorted(glob.glob('solutions/[1-9]*-*.json.submission.result')))
    submission_not_get_result = [x for x in submissions if x + ".result" not in results]

    print(submission_not_get_result)

    for submission in submission_not_get_result:
        with open(submission) as f:
            submission_id = f.read().strip()

        try:
            time.sleep(1.0)
            response = get_submission(submission_id).strip("\"")
            js = json.loads(response)
            failure = "Failure" in js["Success"]["submission"]["score"]
            success = "Success" in js["Success"]["submission"]["score"]
            if failure or success:
                print(f'saving {submission}.result')
                with open(f'{submission}.result', mode="w") as result_file:
                    result_file.write(response)
                    result_file.write("\n")
        except Exception as e:
            print(e)


if __name__ == "__main__":
    main()
