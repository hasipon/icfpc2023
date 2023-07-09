import subprocess

subprocess.run(["cargo", "build", "--release"]);

for n in range(60, 90):
    subprocess.Popen(["./target/release/solver", str(n + 1)])
