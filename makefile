mpi:
	@clear
	mpicc -O3 -I ./ ./polybench.c -o durbin_mpi.out durbin_mpi.c
	@echo Executando programa...
	mpirun -n 6 ./durbin_mpi.out -d small -t 6

pth:
	@clear
	gcc -O3 -I ./ ./polybench.c -o durbin_pthread.out durbin_pthread.c
	@echo Executando programa...
	./durbin_pthread.out -d small -t 6

seq:
	@clear
	gcc -O3 -I ./ ./polybench.c -o durbin.out durbin.c
	@echo Executando programa...
	./durbin.out -d small

clean:
	@clear
	rm -f *.out

run:
	@clear
	@{ time ./durbin.out -d small; }
