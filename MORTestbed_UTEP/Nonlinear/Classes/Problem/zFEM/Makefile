EXE  = zFEM
OBJS = zFEM.o \
       matops.o \
       zfem_elems.o \
       zfem_mass.o \
       zfem_forces.o
CXX = g++44

# C++ Compiler options
CXXFLAGS = -Wall -g

CXXLINKFLAGS = 

all:
	$(CXX) -c matops.cpp zfem_elems.cpp zfem_forces.cpp zfem_mass.cpp
	#$(CXX) $(CXXFLAGS) matops.cpp zfem_elems.cpp zfem_mass.cpp zfem_forces.cpp -lm -o $(EXE)

clean:
	rm matops.o zfem_elems.o zfem_forces.o zfem_mass.o
