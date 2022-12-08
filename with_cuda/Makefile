BIN		:= main

CUDA_PATH 	:= /usr/local/cuda-11.2
INCLUDES 	+= -I. -I$(CUDA_PATH)/include
LIBS		:= -L$(CUDA_PATH)/lib64

CXXFLAGS	:= -O3 -m64 -stdlib=libc++ #-fno-inline -pg
NVCCFLAGS	:= --default-stream per-thread --allow-unsupported-compiler #-g -G #allows per-thread default stream 
LDFLAGS		:= -lgsl -lgslcblas -lgomp 

NVCC		:= nvcc
CXX		:= clang++ #dpcpp
LINKER		:= clang++ #dpcpp

C_SOURCES	:= main.cpp functions.cpp classes.cpp
CU_SOURCES	:= dev_functions.cu
HEADERS		:= functions.hpp classes.hpp dev_functions.hpp

C_OBJS		:= $(patsubst %.cpp, %.o, $(C_SOURCES))
CU_OBJS		:= $(patsubst %.cu, %.o, $(CU_SOURCES))


$(BIN): $(C_OBJS) $(HEADERS) $(CU_OBJS)
	$(LINKER) -o $@ $(C_OBJS) $(LDFLAGS) $(CXXFLAGS) $(LIBS)

$(C_OBJS): $(C_SOURCES) $(HEADERS)
	$(CXX) -c $(C_SOURCES) $(CXXFLAGS)

$(CU_OBJS): $(CU_SOURCES) $(HEADERS)
	$(NVCC) -c $(CU_SOURCES) $(NVCCFLAGS) $(INCLUDES)

clean: 
	rm $(BIN) graphic *.o core*
