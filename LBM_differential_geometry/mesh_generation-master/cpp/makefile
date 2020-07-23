CXX = g++
CXXFLAGS = -I/usr/lib/x86_64-linux-gnu/
INCLUDES = 
LD = icc
LDFLAGS := -L/usr/lib/x86_64-linux-gnu/ -lboost_filesystem -lboost_program_options -lboost_system

#/home/ubuntu/workspace/boost
#/home/ubuntu/workspace/boost/stage/lib

SOURCES = meshgen.cpp
HEADERS = meshgen.h
OBJ     = main.o meshgen.o

.SUFFIXES: .cpp

all: mesh

geometry: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o mesh $(LDFLAGS)
	
main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp
	
diffgeom.o: meshgen.cpp meshgen.h
	$(CXX) $(CXXFLAGS) -c meshgen.cpp

clean:
	rm -f main
	rm -f *.o
	rm -f *.i
	rm -f *.lst
	rm -f output.*