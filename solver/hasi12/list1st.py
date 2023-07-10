import os
import glob
import json
import re
result = {}
for filename in glob.glob("../../solutions/*.json.submission.result"):
	with open(filename) as f:
		data = json.load(f)
		try:
			problem_id = data['Success']['submission']['problem_id']
			score = data['Success']['submission']['score']['Success']
			if problem_id not in result:
				result[problem_id] = (score, filename)
			elif score > result[problem_id][0]:
				result[problem_id] = (score, filename)
		except:
			pass

for x, y in result.items():
	name = y[1][16:-23]
	name2 = re.sub('-hasi12.[0-9]+$', '', name)
	if os.path.exists(f'{name}.err') or os.path.exists(f'{name2}.err'):
		continue
	print(name)
