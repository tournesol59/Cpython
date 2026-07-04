gcc -o libdeviate.o -c libdeviate.c 
ar libdeviate.o libdeviate.a
ranlib libdeviate.a

# dlltool.exe --dllname python3.dll --def python3.def --output-lib libpython.a

