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

Para relizar solo uno de los benchmarks se puede utilizar la funcion
> ./benchmark -t ElectricField


MatrixSub: Se crearán varios archivos en donde se muestra la duración de la operación suma para varios tamaños distintos de matriz.
MatrixAdd: Se crearán varios archivos en donde se muestra la duración de la operación resta para varios tamaños distintos de matriz.
ElectricField: Se generarán 4 gráficas.
	La primera representa el campo vectorial del campo eléctrico generado
	Las otras tres gráficas representan el cambio en magnitud y ángulo de la partícula utilizando un factor de alfa distinto.
------------------------------------------------
Para realizar las pruebas del 
- Doolittle: solver, inverter, substitución hacia adelante y hacia atrás.
- Operaciones simd: suma y resta de matrices
- Malla resistiva: obtener el índice a partir de dos nodos, obtener los nodos a partir del índice de una resistencia y navegar por la malla.
> cd bin
> ./tester

Para realizar únicamente una de las pruebas se puede utilizar el comando 
> ./tester -t ResistorGrid
Esta prueba genera una gráfica con el recorrido de la corriente utilizando las corrientes máximas de cada nodo
------------------------------------------------

