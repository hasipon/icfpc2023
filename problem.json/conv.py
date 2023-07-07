import sys
import json
with open(sys.argv[1]) as f:
	data = json.load(f)
	assert set(data.keys()) == {'Success'}
	print(data['Success'])
