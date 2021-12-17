CXXFLAGS=-O3 -shared -g -Wall -std=c++17 `root-config --cflags` -Iinclude -fPIC
LDFLAGS=`root-config --glibs` -Llib

CC=$(CXX) $(CXXFLAGS) $(LDFLAGS)

all : Directories Utilities LinAlg AggressiveAvocado

Utilities : Directories
	$(CC) -o lib/libUtilities.so src/Utilities.cxx

LinAlg : Directories
	$(CC) -o lib/libLinAlg.so src/LinAlg.cxx

AggressiveAvocado : src/AggressiveAvocado.cxx
	$(CC) -o lib/libAggressiveAvocado.so src/AggressiveAvocado.cxx

Directories :
	mkdir -p lib

clean :
	rm -rf lib/*.so*
