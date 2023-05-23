BIN		:= main

INCLUDES 	+= 
LIBS		:= 

CXXFLAGS	:= -O3 -m64 #-pg
LDFLAGS		:= -std=c++17 -sycl-std=2020 -fsycl -fsycl-targets=nvptx64-nvidia-cuda -Xsycl-target-backend --cuda-gpu-arch=sm_80

CXX		:= clang++ 
LINKER		:= clang++  

C_SOURCES	:= main.cpp functions.cpp classes.cpp
HEADERS		:= functions.hpp classes.hpp

C_OBJS		:= $(patsubst %.cpp, %.o, $(C_SOURCES))


$(BIN): $(C_OBJS) $(HEADERS)
	$(LINKER) -o $@ $(C_OBJS) $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) $(LIBS)

$(C_OBJS): $(C_SOURCES) $(HEADERS)
	$(CXX) -c $(C_SOURCES) $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)

clean: 
	rm $(BIN) graphic *.o core*
