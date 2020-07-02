CXX = g++
CXXFLAGS = -std=c++11 -I/usr/lib/x86_64-linux-gnu/
INCLUDES = 
LD = icc
LDFLAGS := -L/usr/lib/x86_64-linux-gnu/ -lboost_filesystem -lboost_program_options -lboost_system -L/home/ubuntu/workspace/DateTime/Code/

#/home/ubuntu/workspace/boost
#/home/ubuntu/workspace/boost/stage/lib

SOURCES = interp.cpp
HEADERS = interp.h
OBJ     = main.o interp.o

.SUFFIXES: .cpp

all: interpol

interpol: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o interpol $(LDFLAGS)
	
main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp
	
interp.o: interp.cpp interp.h
	$(CXX) $(CXXFLAGS) -c interp.cpp

clean:
	rm -f main
	rm -f *.o
	rm -f *.i
	rm -f *.lst
	rm -f output.*