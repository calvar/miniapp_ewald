BIN		:= main

INCLUDES 	+= 
LIBS		:= 

CXXFLAGS	:= -O3 -m64 -fno-inline #-pg
LDFLAGS		:= -lgsl -lgslcblas -lgomp

CXX		:= g++ #clang++ #dpcpp
LINKER		:= g++ #clang++ #dpcpp

C_SOURCES	:= main.cpp functions.cpp classes.cpp
HEADERS		:= functions.hpp classes.hpp

C_OBJS		:= $(patsubst %.cpp, %.o, $(C_SOURCES))


$(BIN): $(C_OBJS) $(HEADERS)
	$(LINKER) -o $@ $(C_OBJS) $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) $(LIBS) -fopenmp

$(C_OBJS): $(C_SOURCES) $(HEADERS)
	$(CXX) -c $(C_SOURCES) $(CXXFLAGS) $(INCLUDES) -fopenmp

clean: 
	rm $(BIN) graphic *.o core*
