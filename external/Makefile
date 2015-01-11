
all: dirs glowl anttweakbar glew glfw3 optionparser

.PHONY: dirs glowl anttweakbar

dirs:
	@mkdir -p include
	@mkdir -p lib

glowl:
	if cd "glowl"; then echo "already have cloned glowl"; else git clone https://github.com/invor/glowl glowl; fi
	make -C ./glowl/
	cp ./glowl/libglowl.a ./lib/
	cp -r ./glowl/include/glowl ./include/

anttweakbar:
	if cd "anttweakbar-code"; then echo "already have cloned anttweakbar"; else git clone git://git.code.sf.net/p/anttweakbar/code anttweakbar-code; fi
	make -C ./anttweakbar-code/src
	cp ./anttweakbar-code/lib/libAntTweakBar.{a,so*} ./lib/
	cp ./anttweakbar-code/include/AntTweakBar.h ./include

glew:

glfw3:

docopt:
	if cd "docopt.cpp"; then echo "already have cloned docopt"; else git clone https://github.com/docopt/docopt.cpp; fi
	g++ -c ./docopt.cpp/docopt.cpp -o ./docopt.cpp/docopt.o -Wall -g -Wextra -O2 -std=c++11 -fPIC
	ar rcs ./lib/libdocopt.a ./docopt.cpp/docopt.o
	ranlib ./lib/libdocopt.a
	cp ./docopt.cpp/docopt*.h ./include

optionparser:
	wget -P ./include/ http://optionparser.sourceforge.net/optionparser.h
#cp ./optionparser-1.3/src/optionparser.h ./include/

clean:
	rm -rf ./include ./lib

