GET_SPECTRUM_TARGET=GetSpectrum
GET_DESCRIPTOR_TARGET=GetDescriptor
GET_SPECTRUM_SOURCE=GetSpectrum/GetSpectrum.cpp
GET_DESCRIPTOR_SOURCE=GetDescriptor/GetDescriptor.cpp

COMPILER ?= gcc
#COMPILER ?= clang

ifeq ($(COMPILER),gcc)
	CFLAGS += -fopenmp -Wno-deprecated -std=c++14 -pthread -Wno-invalid-offsetof
	LFLAGS += -lgomp -lstdc++ -lpthread
	CC=gcc
	CXX=g++
else
	CFLAGS += -Wno-deprecated -std=c++14 -pthread -Wno-invalid-offsetof -Wno-dangling-else
	LFLAGS += -lstdc++
	CC=clang
	CXX=clang++
endif

CFLAGS += -O3 -DRELEASE -funroll-loops -ffast-math -g
LFLAGS += -O3 -g

BIN = Bin/Linux/
BIN_O = Obj/Linux/
INCLUDE = ./Include/ -I .


MD=mkdir

GET_SPECTRUM_OBJECTS=$(addprefix $(BIN_O), $(addsuffix .o, $(basename $(GET_SPECTRUM_SOURCE))))
GET_DESCRIPTOR_OBJECTS=$(addprefix $(BIN_O), $(addsuffix .o, $(basename $(GET_DESCRIPTOR_SOURCE))))
GET_SPECTRUM_OBJECT_DIR=$(dir $(GET_SPECTRUM_OBJECTS))
GET_DESCRIPTOR_OBJECT_DIR=$(dir $(GET_DESCRIPTOR_OBJECTS))

all: make_dirs
all: $(BIN)$(GET_SPECTRUM_TARGET)
all: $(BIN)$(GET_DESCRIPTOR_TARGET)

getspectrum: make_dirs
getspectrum: $(BIN)$(GET_SPECTRUM_TARGET)

getdescriptor: make_dirs
getdescriptor: $(BIN)$(GET_DESCRIPTOR_TARGET)

clean:
	rm -rf $(BIN)$(GET_SPECTRUM_TARGET)
	rm -rf $(BIN)$(GET_DESCRIPTOR_TARGET)
	rm -rf $(BIN_O)
	cd OpenMesh && make clean
	cd PNG && make clean

make_dirs: FORCE
	$(MD) -p $(BIN)
	$(MD) -p $(BIN_O)
	$(MD) -p $(GET_SPECTRUM_OBJECT_DIR)
	$(MD) -p $(GET_DESCRIPTOR_OBJECT_DIR)

$(BIN)$(GET_SPECTRUM_TARGET): $(GET_SPECTRUM_OBJECTS)
	$(CXX) -o $@ $(GET_SPECTRUM_OBJECTS) -L$(BIN) $(LFLAGS)

$(BIN)$(GET_DESCRIPTOR_TARGET): $(GET_DESCRIPTOR_OBJECTS)
	cd OpenMesh && make COMPILER=$(COMPILER)
	cd PNG && make COMPILER=$(COMPILER)
	$(CXX) -o $@ $(GET_DESCRIPTOR_OBJECTS) -L$(BIN) $(LFLAGS) -lOpenMesh_Core -lecho_PNG -ljpeg -lz

$(BIN_O)%.o: $(SRC)%.cpp
	$(CXX) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<

FORCE: