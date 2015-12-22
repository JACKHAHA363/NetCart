CXX1 = clang++
CXXFLAGS1 = -std=c++1y -stdlib=libc++ -c -g -Wall -Wextra
LD1 = clang++
LDFLAGS1 = -std=c++1y -stdlib=libc++

####################################################
all: bin/netcart_clg bin/netcart_clg_para

###################clang#############################
bin/netcart_clg: bin .objs .objs/main_clg.o .objs/netcart_clg.o .objs/util_clg.o
	$(LD1) .objs/main_clg.o .objs/netcart_clg.o .objs/util_clg.o $(LDFLAGS1) -o bin/netcart_clg

bin/netcart_clg_para: bin .objs .objs/main_clg.o .objs/netcart_clg_para.o .objs/util_clg.o
	$(LD1) .objs/main_clg.o .objs/netcart_clg_para.o .objs/util_clg.o $(LDFLAGS1) -o bin/netcart_clg_para

.objs/main_clg.o: main.cpp
	$(CXX1) $(CXXFLAGS1) main.cpp -o .objs/main_clg.o

.objs/netcart_clg.o: netcart.cpp netcart.h
	$(CXX1) $(CXXFLAGS1) netcart.cpp -o .objs/netcart_clg.o

.objs/netcart_clg_para.o: netcart_para.cpp netcart_para.h
	$(CXX1) $(CXXFLAGS1) netcart_para.cpp -o .objs/netcart_clg_para.o

.objs/util_clg.o: util.cpp util.h
	$(CXX1) $(CXXFLAGS1) util.cpp -o .objs/util_clg.o
############################general##################
bin:
	@mkdir -p bin

.objs:
	@mkdir -p .objs

clean:
	rm -rf bin .objs *.o *.txt
