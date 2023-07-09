import sys
with open(sys.argv[1]) as f:
	a = f.readlines()
	print(len(a))
with open(sys.argv[2], 'w') as f:
	f.write(' '.join(a[-1].split(' ')[1:]))
