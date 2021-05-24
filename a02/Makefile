# -------------------------------------
#  build code

all:
	cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -Wno-dev && cmake --build build --parallel

clean:
	rm -rf build

# -------------------------------------
# build depedency archives

dep: dep_lin dep_win

dep_lin:
	rm -rf depedencies_linux.zip
	zip -r --symlinks depedencies_linux.zip dep/assimp dep/include dep/lin

dep_win:
	rm -rf depedencies_windows.zip
	zip -r depedencies_windows.zip dep/assimp dep/include dep/win
