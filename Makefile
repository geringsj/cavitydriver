
CC=g++
CPPFLAGS=-Wall -g -Wextra -Isrc -std=c++11 
#-Iexternal/include/
LDFLAGS=
#-ldl -lglfw -lGLEW -lGL -lm -lrt -Lexternal/lib/ -Wl,-Rexternal/lib/

SOURCES=$(wildcard src/**/*.cpp src/*.cpp)
OBJECTS=$(patsubst %.cpp,%.o,$(SOURCES))

TESTS_SRC=$(wildcard tests/*.cpp)
TESTS_OBJ=$(patsubst tests/%.cpp,bin/%,$(TESTS_SRC))

MAIN=main.cpp
MAIN_BIN=./bin/main

all: build $(OBJECTS) $(MAIN_BIN)

tests: build $(OBJECTS) $(TESTS_OBJ)

run: all
	$(MAIN_BIN)

opt: CFLAGS += -O3 
#-DNDEBUG
opt: all

build:
	@mkdir -p bin
	@mkdir -p out

$(MAIN_BIN): $(MAIN) $(OBJECTS)
	$(CC) $(CPPFLAGS) -o $(MAIN_BIN) $(MAIN) $(OBJECTS) $(LDFLAGS)

$(TESTS_OBJ): $(TESTS_SRC) $(OBJECTS)
	$(CC) $(CPPFLAGS) -o $@ $(patsubst bin/%,tests/%.cpp,$@) $(filter src%,$^) $(LDFLAGS)

.PHONY: clean
clean:
	rm -rf bin $(OBJECTS)

new: clean all

