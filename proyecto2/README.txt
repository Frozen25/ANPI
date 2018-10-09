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
------------------------------------------------
Se debe contar con un procesador que soporte la tecnología AVX2
------------------------------------------------
Desde la carpeta build se pueden realizar varias tareas:
------------------------------------------------
Para realizar los benchmarks de la resta, para comparar los tiempos con y sin implementar operaciones simd
> cd bin
> ./benchmark
Se crearán varios archivos en donde se muestra la duración de la operación para varios tamaños distintos de matriz.
------------------------------------------------
Para realizar las pruebas del 
- Doolittle: solver, inverter, substitución hacia adelante y hacia atrás.
- Operaciones simd: suma y resta de matrices
- Malla resistiva: obtener el índice a partir de dos nodos, obtener los nodos a partir del índice de una resistencia y navegar por la malla.
> cd bin
> ./tester
------------------------------------------------

