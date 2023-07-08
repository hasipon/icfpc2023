import subprocess

subprocess.run(["cargo", "build", "--release"]);

for n in range(0, 55):
    subprocess.Popen(["./target/release/solver.exe", str(n)])
