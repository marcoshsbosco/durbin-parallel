import os, sys


print("version: seq", file=sys.stderr)
for size in ["small", "medium", "large"]:
    print(f"size: {size}", file=sys.stderr)
    os.system(f"time ./durbin.out -d {size}")

print("version: pth", file=sys.stderr)
for size in ["small", "medium", "large"]:
    print(f"size: {size}", file=sys.stderr)
    for threads in [2, 4, 8]:
        print(f"threads: {threads}", file=sys.stderr)
        os.system(f"time ./durbin_pthread.out -d {size} -t {threads}")

print("version: mpi", file=sys.stderr)
for size in ["small", "medium", "large"]:
    print(f"size: {size}", file=sys.stderr)
    for threads in [2, 4, 8]:
        print(f"threads: {threads}", file=sys.stderr)
        os.system(f"time mpirun --oversubscribe -n {threads} ./durbin_mpi.out -d {size}")

print("version: hyb", file=sys.stderr)
for size in ["small", "medium", "large"]:
    print(f"size: {size}", file=sys.stderr)
    for threads in [(2, 1), (2, 2), (2, 4)]:
        print(f"threads: {threads}", file=sys.stderr)
        os.system(f"time mpirun --oversubscribe -n {threads[0]} ./durbin_mpi.out -d {size} -t {threads[1]}")
