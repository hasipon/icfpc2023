import os
import glob
import json
import re
result = {}
maxscore = {}
for filename in glob.glob("../../solutions/*.json.submission.result"):
	with open(filename) as f:
		data = json.load(f)
		try:
			problem_id = data['Success']['submission']['problem_id']
			score = data['Success']['submission']['score']['Success']
			if problem_id not in result:
				result[problem_id] = (score, filename)
				maxscore[problem_id] = 0
			elif score > result[problem_id][0]:
				result[problem_id] = (score, filename)
		except:
			pass



for filename in glob.glob("*.err"):
	with open(filename) as f:
		for line in list(f)[::-1]:
			if line.startswith('score = '):
				problem_id = int(filename.split('-')[0])
				maxscore[problem_id] = max(maxscore[problem_id], int(line[8:]))
				break

output = []
for x, y in result.items():
	output.append((maxscore[x] - y[0], x))

for x in sorted(output):
	print(x)
