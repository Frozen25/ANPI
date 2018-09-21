# ANPI-Tarea4

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

Para correr el benchmark para medir tiempos se debe correr
el test benchmarkMatrixAdd.cpp ejecutando 
BOOST_AUTO_TEST_SUITE( MatrixDescomposition )

Los archivos de texto generados por la prueba de en tiempo
se guardan en cmake-build-debug/bin.

Para correr los test LU se debe ingresar al documento LUtest.cpp
ejecuntando el test BOOST_AUTO_TEST_SUITE( LU ).

En el archivo solve.hpp se encuentra los metodos invert y los metodos
solve para los algoritmos LU.

