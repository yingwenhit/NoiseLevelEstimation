call vcvars32
mkdir MakefilesMSVC
cd    MakefilesMSVC
cmake -G "NMake Makefiles" -D EXECUTABLE_OUTPUT_PATH:PATH=".." ..
nmake -f Makefile
cd ..
