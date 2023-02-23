BIN		:= main

CUDA_PATH 	:= #/usr/local/cuda-11.2
INCLUDES 	+= #-I/usr/local/cuda/include
LIBS		:= #-L/usr/local/cuda/lib64

CXXFLAGS	:= -O3 -m64 #-stdlib=libc++ #-fno-inline -pg
NVCCFLAGS	:= -O2 #--default-stream per-thread #--allow-unsupported-compiler #-g -G #allows per-thread default stream 
LDFLAGS		:= #-lgomp #-lgsl -lgslcblas 

NVCC		:= nvcc
CXX		:= g++ #clang++ #dpcpp
LINKER		:= nvcc

C_SOURCES	:= functions.cpp classes.cpp
CU_SOURCES	:= dev_functions.cu
HEADERS		:= functions.hpp classes.hpp dev_functions.hpp 

C_OBJS		:= $(patsubst %.cpp, %.o, $(C_SOURCES))
CU_OBJS		:= $(patsubst %.cu, %.o, $(CU_SOURCES))


$(BIN): $(C_OBJS) $(CU_OBJS) $(HEADERS) 
	$(LINKER) -o $@ main.cpp $(C_OBJS) $(CU_OBJS) $(LDFLAGS) $(NVCCFLAGS) $(LIBS)

$(C_OBJS): $(C_SOURCES) $(HEADERS)
	$(CXX) -c $(C_SOURCES) $(CXXFLAGS)

$(CU_OBJS): $(CU_SOURCES) $(HEADERS)
	$(NVCC) -c $(CU_SOURCES) $(NVCCFLAGS) $(INCLUDES)


clean: 
	rm $(BIN) *.o core*
