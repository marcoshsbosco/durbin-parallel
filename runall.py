import os


print("version: seq")
for size in ["small", "medium", "large"]:
    print(f"size: {size}")
    os.system(f"time ./durbin.out -d {size}")

print("version: pth")
for size in ["small", "medium", "large"]:
    print(f"size: {size}")
    for threads in [2, 4, 8]:
        print(f"threads: {threads}")
        os.system(f"time ./durbin_pthread.out -d {size} -t {threads}")

print("version: mpi")
for size in ["small", "medium", "large"]:
    print(f"size: {size}")
    for threads in [2, 4, 8]:
        print(f"threads: {threads}")
        os.system(f"time mpirun --oversubscribe -n {threads} ./durbin_mpi.out -d {size}")

print("version: hyb")
for size in ["small", "medium", "large"]:
    print(f"size: {size}")
    for threads in [(2, 1), (2, 2), (2, 4)]:
        print(f"threads: {threads}")
        os.system(f"time mpirun --oversubscribe -n {threads[0]} ./durbin_mpi.out -d {size} -t {threads[1]}")
