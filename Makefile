
CXX=g++
#CXX=clang++
CPPFLAGS=-Wall -g -Wextra -Isrc -std=c++11 -Iexternal/include/
LDFLAGS= -ldl -lglowl -lAntTweakBar -lglfw -lGLEW -lGL -lm -lrt -Lexternal/lib/ -Wl,-Rexternal/lib/

SOURCES=$(wildcard src/**/*.cpp src/*.cpp)
OBJECTS=$(patsubst %.cpp,%.o,$(SOURCES))

TESTS_SRC=$(wildcard tests/*.cpp)
TESTS_OBJ=$(patsubst tests/%.cpp,bin/%,$(TESTS_SRC))

MAINS=$(wildcard cavity*.cpp)
MAINS_BIN=$(patsubst %.cpp,bin/%,$(MAINS))

all: build $(OBJECTS) $(MAINS_BIN)

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

build:
	@mkdir -p bin
	@mkdir -p out
#@mkdir -p build

$(MAINS_BIN): $(MAINS) $(OBJECTS)
	$(CXX) $(CPPFLAGS) -o $@ $(patsubst bin/%,%.cpp,$@) $(OBJECTS) $(LDFLAGS)

$(TESTS_OBJ): $(TESTS_SRC) $(OBJECTS)
	$(CXX) $(CPPFLAGS) -o $@ $(patsubst bin/%,tests/%.cpp,$@) $(filter src%,$^) $(LDFLAGS)

doxy:
	doxygen ./doxygenconfig.txt

.PHONY: clean
clean:
	rm -rf bin $(OBJECTS) src/*.hpp.gch documentation 
#build

new: clean all

