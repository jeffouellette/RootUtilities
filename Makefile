CXXFLAGS=-O3 -shared -g -Wall `root-config --cflags` -Iinclude -fPIC
LDFLAGS=`root-config --glibs` -Llib

CC=$(CXX) $(CXXFLAGS) $(LDFLAGS)

all : Utilities

Utilities :
	$(CC) -o lib/libUtilities.so src/Utilities.cxx

clean :
	rm -rf lib/*.so*
