TARGET = echo_PNG
SOURCE = \
	png.c \
	pngerror.c \
	pnggccrd.c \
	pngget.c \
	pngmem.c \
	pngpread.c \
	pngread.c \
	pngrio.c \
	pngrtran.c \
	pngrutil.c \
	pngset.c \
	pngtest.c \
	pngtrans.c \
	pngvcrd.c \
	pngwio.c \
	pngwrite.c \
	pngwtran.c \
	pngwutil.c
	
COMPILER ?= gcc

TARGET_LIB = lib$(TARGET).a

ifeq ($(COMPILER),gcc)
	CFLAGS += -fopenmp -Wno-deprecated -Wno-unused-result
	CC=gcc
	CXX=g++
else
	CFLAGS += -Wno-deprecated -Wno-dangling-else -Wno-unused-result
	CC=clang
	CXX=clang++
endif

CFLAGS += -O3 -DRELEASE -funroll-loops -ffast-math -DNDEBUG

BIN = ../Bin/Linux/
BIN_O = ../Obj/Linux/$(TARGET)/
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

$(BIN)$(TARGET_LIB): $(OBJECTS)
	$(AR) rcs $@ $(OBJECTS)

$(BIN_O)%.o: %.c
	$(CC) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<

FORCE:
