
all: dirs glm glowl glew glfw3 anttweakbar optionparser

.PHONY: dirs glowl anttweakbar

dirs:
	@mkdir -p include
	@mkdir -p lib

glowl: dirs
	if cd "glowl"; then echo "already have cloned glowl"; else git clone https://github.com/invor/glowl glowl; fi
	cp -r ./glm/glm ./glowl/include/glowl
	make -C ./glowl/
	cp ./glowl/libglowl.a ./lib/
	cp -r ./glowl/include/glowl ./include/

anttweakbar: dirs
	if cd "anttweakbar-code"; then echo "already have cloned anttweakbar"; else git clone git://git.code.sf.net/p/anttweakbar/code anttweakbar-code; fi
	make -C ./anttweakbar-code/src
	cp ./anttweakbar-code/lib/libAntTweakBar.{a,so*} ./lib/
	cp ./anttweakbar-code/include/AntTweakBar.h ./include

glew: dirs
	if cd "glew"; then echo "already have cloned glew"; else git clone https://github.com/nigels-com/glew.git glew; fi
	make -C glew extensions
	make -C glew glew.lib
	cp glew/lib/libGLEW.so* ./lib
	cp -r glew/include/GL ./include

glfw3: dirs
	if cd "glfw"; then echo "already have cloned glfw"; else git clone https://github.com/glfw/glfw glfw; fi
	cd glfw && cmake -DBUILD_SHARED_LIBS=ON
	make -C glfw sharing
	cp ./glfw/src/libglfw* ./lib
	cp -r ./glfw/include/GLFW ./include

glm: dirs
	if cd "glm"; then echo "already have cloned glm"; else git clone https://github.com/g-truc/glm glm; fi
	cp -r ./glm/glm ./include

optionparser: dirs
	wget -P ./include/ http://optionparser.sourceforge.net/optionparser.h
#cp ./optionparser-1.3/src/optionparser.h ./include/

clean:
	rm -rf ./include ./lib

