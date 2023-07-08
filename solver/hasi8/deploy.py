import os
import glob
import json
for filename in glob.glob("*.out"):
	with open(filename) as f:
		a = list(f)
		for i in range(0, len(a)):
			try:
				json.loads(a[-(i+1)])
				new_path = f'../../solutions/{filename[:-4]}-hasi8.{len(a)}.json'
				if not os.path.exists(new_path):
					with open(new_path, 'w') as out:
						out.write(a[-(i+1)])
				break
			except:
				pass
