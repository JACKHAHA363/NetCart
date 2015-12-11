CXX1 = clang++
CXXFLAGS1 = -std=c++1y -stdlib=libc++ -c -g -Wall -Wextra
LD1 = clang++
LDFLAGS1 = -std=c++1y -stdlib=libc++

CXX2 = g++
CXXFLAGS2 = -std=c++11 -c -g -Wall -Wextra
LD2 = g++
LDFLAGS2 = -std=c++11

all: netcart_clg netcart_gcc

###################clang#############################
netcart_clg: bin .objs main_clg.o netcart_clg.o util_clg.o
	$(LD1) .objs/main_clg.o .objs/netcart_clg.o .objs/util_clg.o $(LDFLAGS1) -o bin/netcart_clg

main_clg.o: main.cpp
	$(CXX1) $(CXXFLAGS1) main.cpp -o .objs/main_clg.o

netcart_clg.o: netcart.cpp netcart.h
	$(CXX1) $(CXXFLAGS1) netcart.cpp -o .objs/netcart_clg.o

util_clg.o: util.cpp util.h
	$(CXX1) $(CXXFLAGS1) util.cpp -o .objs/util_clg.o

##################gcc###############################
netcart_gcc: bin .objs main_gcc.o netcart_gcc.o util_gcc.o
	$(LD2) .objs/main_gcc.o .objs/netcart_gcc.o .objs/util_gcc.o $(LDFLAGS2) -o bin/netcart_gcc

main_gcc.o: main.cpp
	$(CXX2) $(CXXFLAGS2) main.cpp -o .objs/main_gcc.o

netcart_gcc.o: netcart.cpp netcart.h
	$(CXX2) $(CXXFLAGS2) netcart.cpp -o .objs/netcart_gcc.o

util_gcc.o: util.cpp util.h
	$(CXX2) $(CXXFLAGS2) util.cpp -o .objs/util_gcc.o

############################general##################
bin:
	mkdir bin

.objs:
	mkdir .objs

clean:
	rm -rf bin .objs *.o *.txt
