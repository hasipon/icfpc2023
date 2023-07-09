import os
import glob
import json
import urllib.request
import urllib.parse
import urllib.error
import time


def is_valid_json(s):
    try:
        json.loads(s)
    except ValueError as err:
        return False
    return True


def main():
    for solution in glob.glob('solutions/[1-9]*-*.json'):
        with open(solution) as f:
            s = f.read()
            if not is_valid_json(s):
                print("remove invalid json:", solution, "content:", s)
                os.remove(solution)


if __name__ == "__main__":
    main()
