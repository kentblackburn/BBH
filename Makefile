# Makefile for the Binary Black Hole Code
# James Kent Blackburn
# Oct 18th, 1993

CFLAGS = -xO4 

CFILES = bracket.c brent.c evaluate.c fastvol.c scalar.c sort.c

newext: newextreme.c $(CFILES)
	cc -o newext $(CFLAGS) $(CFILES) newextreme.c -lm

clean: 
	rm -f *.o *.u
