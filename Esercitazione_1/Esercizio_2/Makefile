CC = g++
CFLAGS = -Wall -O3 --std=c++11

Esercizio_1.exe : main.o random.o
	$(CC) random.o main.o -o main.exe
Esercizio_1.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Esercizio_1.exe seed.out
