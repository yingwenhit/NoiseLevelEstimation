mkdir MakefilesMinGW
cd    MakefilesMinGW
cmake -G "MinGW Makefiles" -D EXECUTABLE_OUTPUT_PATH:PATH=".." ..
mingw32-make -f Makefile
cd ..
