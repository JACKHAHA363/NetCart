CXX = clang++
CXXFLAGS = -std=c++1y -stdlib=libc++ -c -g -Wall -Wextra
LD = clang++
LDFLAGS = -std=c++1y -stdlib=libc++

netcart: main.o netcart.o util.o
	$(LD) main.o netcart.o util.o $(LDFLAGS) -o netcart

main.o: main.cpp
	$(CXX) $(CXXFLAGS) main.cpp

netcart.o: netcart.cpp netcart.h
	$(CXX) $(CXXFLAGS) netcart.cpp

util.o: util.cpp util.h
	$(CXX) $(CXXFLAGS) util.cpp

clean:
	rm -rf *.o netcart