make:
	@clear
	gcc -O3 -I ./ ./polybench.c -o durbin_pthread.out durbin_pthread.c
	@echo Executando programa...
	./durbin_pthread.out -d small -t 1

clean:
	@clear
	rm -f *.out

run:
	@clear
	@{ time ./durbin.out -d small; }

seq:
	@clear
	gcc -O3 -I ./ ./polybench.c -o durbin.out durbin.c
	@echo Executando programa...
	./durbin.out -d small
