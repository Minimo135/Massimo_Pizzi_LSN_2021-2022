CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11 -g

Esercizio_10.exe : main.o random.o Genetics.o global.o
	$(CC) random.o main.o Genetics.o global.o -o Esercizio_10.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
global.o: global.h global.cpp
	$(CC) -c global.cpp -o global.o $(CFLAGS)
Genetics.o: Genetics.h Genetics.cpp random.h random.cpp global.h
	$(CC) -c Genetics.cpp -o Genetics.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Esercizio_9.exe seed.out
