CC = g++
CFLAGS = -Wall -O3 --std=c++17

main.exe : main.o random.o stats.o
	$(CC) random.o stats.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : ../random.cpp ../random.h
	$(CC) -c ../random.cpp -o random.o $(CFLAGS)
stats.o : ../stats.cpp ../stats.h
	$(CC) -c ../stats.cpp -o stats.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
