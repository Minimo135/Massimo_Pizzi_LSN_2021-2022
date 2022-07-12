CC = g++
CFLAGS = -Wall -O3 -g

Esercizio_9_ALT.exe : main_ALT.o random.o Genetics.o global.o
	$(CC) random.o main_ALT.o Genetics.o global.o -o Esercizio_9_ALT.exe
main_ALT.o : main_ALT.cpp
	$(CC) -c main_ALT.cpp -o main_ALT.o $(CFLAGS)

Esercizio_9.exe : main.o random.o Genetics.o global.o
	$(CC) random.o main.o Genetics.o global.o -o Esercizio_9.exe
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
erase:
	sh clean.sh
