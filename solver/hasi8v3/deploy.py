import os
import glob
import json
for filename in glob.glob("[0-9]*.out"):
	with open(filename) as f:
		a = list(f)
		for i in range(0, len(a)):
			try:
				json.loads(a[-(i+1)])
				new_path = f'../../solutions/{filename[:-4]}-hasi8v3.{len(a)-i}.json'
				if not os.path.exists(new_path):
					with open(new_path, 'w') as out:
						out.write(a[-(i+1)])
				break
			except:
				pass
