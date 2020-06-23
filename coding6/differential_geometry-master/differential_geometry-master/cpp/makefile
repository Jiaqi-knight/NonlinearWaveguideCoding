CXX = g++
CXXFLAGS = -I/usr/lib/x86_64-linux-gnu/ -I/home/ubuntu/workspace/Finite_differences/Code -I/home/ubuntu/workspace/Determinant/Code
INCLUDES = 
LD = icc
LDFLAGS := -L/usr/lib/x86_64-linux-gnu/ -lboost_filesystem -lboost_program_options -lboost_system

#/home/ubuntu/workspace/boost
#/home/ubuntu/workspace/boost/stage/lib

SOURCES = diffgeom.cpp
HEADERS = diffgeom.h
OBJ     = main.o diffgeom.o

.SUFFIXES: .cpp

all: geometry

geometry: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o geometry $(LDFLAGS)
	
main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp
	
diffgeom.o: diffgeom.cpp diffgeom.h
	$(CXX) $(CXXFLAGS) -c diffgeom.cpp

clean:
	rm -f main
	rm -f *.o
	rm -f *.i
	rm -f *.lst
	rm -f output.*