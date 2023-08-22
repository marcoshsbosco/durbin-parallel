make:
	@clear
	gcc -O3 -I ./ ./polybench.c -o durbin.out durbin.c
	@echo Executando programa...
	@{ time ./durbin.out -h; }

clean:
	@clear
	rm -f *.out
