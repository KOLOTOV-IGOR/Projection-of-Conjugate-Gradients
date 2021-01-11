CC = gcc
CFLAGS = -g -Wall -m64
LDLIBS = -lm

mem.o: memory.c memory.h
	$(CC) $(CFLAGS) -c memory.c -o mem.o

aux.o: auxelirable.c auxelirable.h
	$(CC) $(CFLAGS) -c auxelirable.c $(LDLIBS) -o aux.o

in_eq.o: IntegralEquation/IntegralEquation.c IntegralEquation/IntegralEquation.h
	$(CC) $(CFLAGS) -c IntegralEquation/IntegralEquation.c $(LDLIBS) -o in_eq.o 

ker.o: IntegralEquation/kernel_funcs.c IntegralEquation/kernel_funcs.h
	$(CC) $(CFLAGS) -c IntegralEquation/kernel_funcs.c $(LDLIBS) -o ker.o 

prog: main.c in_eq.o ker.o aux.o mem.o
	$(CC) $(CFLAGS) main.c mem.o aux.o ker.o in_eq.o $(LDLIBS) -o prog

run: prog
	./prog && python3 graphics.py

clean:
	rm -f *.o prog
