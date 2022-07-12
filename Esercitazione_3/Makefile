CC = g++
CFLAGS = -Wall -O3 --std=c++11

Esercizio_3.exe : main.o random.o
	$(CC) random.o main.o -o Esercizio_3.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o *.dat Esercizio_3.exe seed.out