CXX = clang++
CXXFLAGS = -std=c++1y -stdlib=libc++ -c -g -Wall -Wextra
LD = clang++
LDFLAGS = -std=c++1y -stdlib=libc++

netcart: bin .objs main.o netcart.o util.o
	$(LD) .objs/main.o .objs/netcart.o .objs/util.o $(LDFLAGS) -o bin/netcart

main.o: main.cpp
	$(CXX) $(CXXFLAGS) main.cpp -o .objs/main.o

netcart.o: netcart.cpp netcart.h
	$(CXX) $(CXXFLAGS) netcart.cpp -o .objs/netcart.o

util.o: util.cpp util.h
	$(CXX) $(CXXFLAGS) util.cpp -o .objs/util.o

bin:
	mkdir bin

.objs:
	mkdir .objs

clean:
	rm -rf bin .objs *.o *.txt
