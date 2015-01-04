
CXX=g++
#CXX=clang++
CPPFLAGS=-Wall -g -Wextra -Isrc -std=c++11 -Iexternal/include/
LDFLAGS=
#-ldl -lglowl -lAntTweakBar -lglfw -lGLEW -lGL -lm -lrt -Lexternal/lib/ -Wl,-Rexternal/lib/

SOURCES=$(wildcard src/**/*.cpp src/*.cpp)
OBJECTS=$(patsubst %.cpp,%.o,$(SOURCES))

TESTS_SRC=$(wildcard tests/*.cpp)
TESTS_OBJ=$(patsubst tests/%.cpp,bin/%,$(TESTS_SRC))

DRIVER=cavitydriver.cpp
DRIVER_BIN=bin/cavitydriver
DRIVER_OBJECTS=$(filter-out src/%Renderer.o src/%System.o,$(OBJECTS)) 

BAKER=cavitybaker.cpp
BAKER_BIN=bin/cavitybaker
BAKER_OBJECTS=$(filter src/%eters.o src/Bakery.o src/Boundary.o src/GridFunction.o src/%Renderer.o src/%System.o src/%eters.o,$(OBJECTS)) 

PAINTER=cavitypainter.cpp
PAINTER_BIN=bin/cavitypainter
PAINTER_OBJECTS=$(filter src/%Renderer.o src/%System.o src/%eters.o,$(OBJECTS)) 

.PHONY: clean painter

all: build $(DRIVER_BIN)

tests: build $(OBJECTS) $(TESTS_OBJ)

run: all
	$(MAINS_BIN)

mpi: CXX = mpic++
mpi: CPPFLAGS += -DWITHMPI
mpi: all

optmpi: CXX = mpic++
optmpi: CPPFLAGS += -DWITHMPI -O3 -flto -fwhole-program -DNDEBUG
optmpi: all

opt: CPPFLAGS += -O3 -flto -fwhole-program -DNDEBUG
opt: all

painter: LDFLAGS += -ldl -lglowl -lAntTweakBar -lglfw -lGLEW -lGL -lm -lrt -Lexternal/lib/ -Wl,-Rexternal/lib/
painter: build $(PAINTER_BIN)

baker: LDFLAGS += -ldl -ldocopt -lglowl -lAntTweakBar -lglfw -lGLEW -lGL -lm -lrt -Lexternal/lib/ -Wl,-Rexternal/lib/
baker: build $(BAKER_BIN)

build:
	@mkdir -p bin
	@mkdir -p out
#@mkdir -p build

$(DRIVER_BIN): $(DRIVER) $(DRIVER_OBJECTS)
	$(CXX) $(CPPFLAGS) -o $(DRIVER_BIN) $(DRIVER) $(DRIVER_OBJECTS) $(LDFLAGS)

$(PAINTER_BIN): $(PAINTER) $(PAINTER_OBJECTS)
	$(CXX) $(CPPFLAGS) -o $(PAINTER_BIN) $(PAINTER) $(PAINTER_OBJECTS) $(LDFLAGS)

$(BAKER_BIN): $(BAKER) $(BAKER_OBJECTS)
	$(CXX) $(CPPFLAGS) -o $(BAKER_BIN) $(BAKER) $(BAKER_OBJECTS) $(LDFLAGS)

$(TESTS_OBJ): $(TESTS_SRC) $(OBJECTS)
	$(CXX) $(CPPFLAGS) -o $@ $(patsubst bin/%,tests/%.cpp,$@) $(filter src%,$^) $(LDFLAGS)

doxy:
	doxygen ./doxygenconfig.txt

clean:
	rm -rf bin $(OBJECTS) src/*.hpp.gch documentation 
#build

new: clean all

