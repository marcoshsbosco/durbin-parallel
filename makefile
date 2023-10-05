all:
	@clear
	gcc -O3 -I ./ ./polybench.c -o durbin.out durbin.c
	gcc -O3 -I ./ ./polybench.c -o durbin_pthread.out durbin_pthread.c
	mpicc -O3 -I ./ ./polybench.c -o durbin_mpi.out durbin_mpi.c
	mpicc -O3 -I ./ ./polybench.c -o durbin_pthread_mpi.out durbin_pthread_mpi.c

hyb:
	@clear
	mpicc -O3 -I ./ ./polybench.c -o durbin_pthread_mpi.out durbin_pthread_mpi.c

mpi:
	@clear
	mpicc -O3 -I ./ ./polybench.c -o durbin_mpi.out durbin_mpi.c

pth:
	@clear
	gcc -O3 -I ./ ./polybench.c -o durbin_pthread.out durbin_pthread.c

seq:
	@clear
	gcc -O3 -I ./ ./polybench.c -o durbin.out durbin.c

clean:
	@clear
	rm -f *.out
