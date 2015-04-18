CC=g++
CFLAGS=-O4

hh: HodgkinHuxley.o RungeKutta.o
	$(CC) -o rel -D HH reliability.cpp HodgkinHuxley.o RungeKutta.o

ml: MorrisLecar.o RungeKutta.o
	$(CC) -o rel -D ML reliability.cpp MorrisLecar.o RungeKutta.o

clean:
	rm *.o
