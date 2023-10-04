import os


print("version: seq")
for size in ["small", "medium", "large"]:
    print(f"size: {size}")
    os.system(f"time ./durbin.out -d {size}")

print("version: pth")
for size in ["small", "medium", "large"]:
    print(f"size: {size}")
    for threads in range(2, 9, 2):
        print(f"threads: {threads}")
        os.system(f"time ./durbin_pthread.out -d {size} -t {threads}")

print("version: mpi")
for size in ["small", "medium", "large"]:
    print(f"size: {size}")
    for threads in range(2, 9, 2):
        print(f"threads: {threads}")
        os.system(f"time mpirun --oversubscribe -n {threads} ./durbin_mpi.out -d {size}")
