CC = g++
CFLAGS = -Wall -O3 --std=c++11

Esercizio_4.exe : MD_MC.o random.o
	$(CC) random.o MD_MC.o -o Esercizio_4.exe
MD_MC.o : MD_MC.cpp
	$(CC) -c MD_MC.cpp -o MD_MC.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Esercizio_4.exe
	sh clean.sh
