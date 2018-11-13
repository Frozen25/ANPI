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
Para utilizar la interfaz mediante consola, se utilizan los flags asignados en la documentación (también se puede ver usando --help)
> cd bin
> ./proyecto3 --help

También, se cuentan con archivos con perfiles de temperatura predefinidos que puede ser llamados de la forma

> ./proyecto3 -p ../../data/perfil.txt

Esto sirve para obtener los valores de temperatura en la placa. Para observar las gráficas de ellas, se corre únicamente el siguiente comandos puesto que las demás banderas definidas
con anterioridad han sido guardadas para este proceso.

> ./benchmark -t thermalGrid

------------------------------------------------

