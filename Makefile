CXX1 = g++
CXXFLAGS1 = -std=c++1y -c -g -Wall -Wextra
LD1 = g++
LDFLAGS1 = -std=c++1y -pthread

####################################################
all: bin/netcart

####################################################
bin/netcart: bin .objs .objs/main.o .objs/netcart.o .objs/util.o
	$(LD1) .objs/main.o .objs/netcart.o .objs/util.o $(LDFLAGS1) -o bin/netcart

.objs/main.o: main.cpp
	$(CXX1) $(CXXFLAGS1) main.cpp -o .objs/main.o

.objs/netcart.o: netcart.cpp netcart.h
	$(CXX1) $(CXXFLAGS1) netcart.cpp -o .objs/netcart.o

.objs/util.o: util.cpp util.h
	$(CXX1) $(CXXFLAGS1) util.cpp -o .objs/util.o

#####################################################
bin:
	@mkdir -p bin

.objs:
	@mkdir -p .objs

clean:
	rm -rf bin .objs *.o *.txt
