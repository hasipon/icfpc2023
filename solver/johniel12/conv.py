import sys
import json
with open(sys.argv[1]) as f:
	data = json.load(f)
	print(len(data['placements']))
	for v in data['placements']:
		print(v['x'], v['y'])
