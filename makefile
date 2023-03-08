CC=g++

CPPFLAGS = -c -g -Wall -O3 -std=c++17
HEADER_PATH = -I "./seqan/include" -I "./boost"

alphaHiASM: main.o overlap.o assembly.o utility.o
	$(CC) -o $@ $^ -lpthread
	
main.o: main.cpp
	$(CC) $(CPPFLAGS) main.cpp $(HEADER_PATH)

overlap.o: overlap.cpp 
	$(CC) $(CPPFLAGS) overlap.cpp $(HEADER_PATH)

assembly.o: assembly.cpp
	$(CC) $(CPPFLAGS) assembly.cpp $(HEADER_PATH)
	
utility.o: utility.cpp
	$(CC) $(CPPFLAGS) utility.cpp $(HEADER_PATH)


all: alphaHiASM
	rm -f *.o

clean:
	rm -f *.o
	rm alphaHiASM