.PHONY: clean

fplot: fplot.c config.h functions.h
	cc -o $@ fplot.c -lnotcurses-core -lm

clean:
	rm -f fplot
