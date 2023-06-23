BIN		:= main

INCLUDES 	+= 
LIBS		:= 

CXXFLAGS	:= -O3 -m64 #-pg
LDFLAGS		:= --stdlib=libc++ -fopenmp -fopenmp-targets=nvptx64 -Xopenmp-target=nvptx64 -march=sm_80 -Rpass=openmp-opt -Rpass-missed=openmp-opt -Rpass-analysis=openmp-opt -fopenmp-cuda-mode -fopenmp-assume-no-nested-parallelism -fopenmp-assume-no-thread-state#-foffload="-lm" -lm -fopenmp -fcf-protection=none -fno-stack-protector -no-pie  #-lgsl -lgslcblas 

CXX		:= clang++ #g++
LINKER		:= clang++ #g++ 

C_SOURCES	:= main.cpp functions.cpp classes.cpp
HEADERS		:= functions.hpp classes.hpp

C_OBJS		:= $(patsubst %.cpp, %.o, $(C_SOURCES))


$(BIN): $(C_OBJS) $(HEADERS)
	$(LINKER) -o $@ $(C_OBJS) $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) $(LIBS)

$(C_OBJS): $(C_SOURCES) $(HEADERS)
	$(CXX) -c $(C_SOURCES) $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)

clean: 
	rm $(BIN) graphic *.o core*
