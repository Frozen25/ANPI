Create a directory build:

> mkdir build;

Go into that directory

> cd build;

You can choose to build a release version with:

> cmake ../ -DCMAKE_BUILD_TYPE=Release

or a debug version with

> cmake ../ -DCMAKE_BUILD_TYPE=Debug

And build everything with

> make

The executables will be stored at build/bin.


To execute the benchmarks you will need python2.7 and python-tk

> sudo apt install python2.7 python-tk

Additionally, you need matplotlib in python2.7

> pip install --user matplotlib

Para ejecutar deflación, se pasa por argumento un polinomio usando bmt::polynomial<T> poly; poly = {c, b a},
donde c, b, a son los coeficientes del polinomio y T un tipo cualquiera, el polinomio se pasa en el parámetro
poly, la raíz en root y debe ser tipo T y el residuo del tipo del polinomio.
