make:
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
