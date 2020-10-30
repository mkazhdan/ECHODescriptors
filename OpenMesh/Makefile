TARGET = OpenMesh_Core
SOURCE = \
	Core/IO/BinaryHelper.cc \
	Core/IO/IOManager.cc \
	Core/IO/OMFormat.cc \
	Core/IO/reader/BaseReader.cc \
	Core/IO/reader/OBJReader.cc \
	Core/IO/reader/OFFReader.cc \
	Core/IO/reader/OMReader.cc \
	Core/IO/reader/PLYReader.cc \
	Core/IO/reader/STLReader.cc \
	Core/IO/writer/BaseWriter.cc \
	Core/IO/writer/OBJWriter.cc \
	Core/IO/writer/OFFWriter.cc \
	Core/IO/writer/OMWriter.cc \
	Core/IO/writer/PLYWriter.cc \
	Core/IO/writer/STLWriter.cc \
	Core/Mesh/ArrayKernel.cc \
	Core/Mesh/ArrayKernelT.cc \
	Core/Mesh/BaseKernel.cc \
	Core/Mesh/PolyConnectivity.cc \
	Core/Mesh/TriConnectivity.cc \
	Core/System/omstream.cc \
	Core/Utils/BaseProperty.cc \
	Core/Utils/Endian.cc
	
COMPILER ?= gcc

TARGET_LIB = lib$(TARGET).a

ifeq ($(COMPILER),gcc)
	CFLAGS += -fopenmp -Wno-deprecated -std=c++14 -pthread -Wno-invalid-offsetof -Wno-unused-result
	CC=gcc
	CXX=g++
else
	CFLAGS += -Wno-deprecated -std=c++14 -pthread -Wno-invalid-offsetof -Wno-dangling-else -Wno-unused-result
	CC=clang
	CXX=clang++
endif

CFLAGS += -O3 -DRELEASE -funroll-loops -ffast-math -DNDEBUG

BIN = ../Bin/Linux/
BIN_O = ../Obj/Linux/
INCLUDE = ../

MD = mkdir
AR = ar

OBJECTS=$(addprefix $(BIN_O), $(addsuffix .o, $(basename $(SOURCE))))
OBJECT_DIRS=$(dir $(OBJECTS))

all: make_dirs
all: $(BIN)$(TARGET_LIB)

clean:
	rm -rf $(BIN_O)
	rm -f $(BIN)$(TARGET_LIB)

make_dirs: FORCE
	$(MD) -p $(BIN)
	$(MD) -p $(BIN_O)
	for dir in $(OBJECT_DIRS); do $(MD) -p $$dir; done

$(BIN)$(TARGET_LIB): $(OBJECTS)
	$(AR) rcs $@ $(OBJECTS)

$(BIN_O)%.o: %.cc
	$(CXX) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<

FORCE:
