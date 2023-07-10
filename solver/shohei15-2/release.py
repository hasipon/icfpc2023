import subprocess

subprocess.run(["cargo", "build", "--release"]);

for n in range(16, 23):
    subprocess.Popen(["./target/release/solver", str(n + 1)])
