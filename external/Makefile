
all: dirs glowl anttweakbar glew glfw3

.PHONY: dirs glowl anttweakbar

dirs:
	@mkdir -p include
	@mkdir -p lib

glowl:
	if cd "glowl"; then echo "already have cloned glowl"; else git clone git://zuggerwad.de:~/glowl glowl; fi
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

clean:
	rm -rf ./include ./lib
