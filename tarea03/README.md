# Tarea #3 - ANPI

#Tarea #3 De Análisis Numérico para Ingeniería
#*Recomendable usar CLion y Linux para manejar el código

#Se puede instar Boost usando

> sudo apt install libboost-all-dev

#Si el boost aún así no funciona, pueden usar:

> sudo apt-get install cmake libblkid-dev e2fslibs-dev libboost-all-dev libaudit-dev

#Create a directory build:

> mkdir build;

#Go into that directory

> cd build;

#You can choose to build a release version with:

> cmake ../ -DCMAKE_BUILD_TYPE=Release

#or a debug version with

> cmake ../ -DCMAKE_BUILD_TYPE=Debug

#And build everything with

> make

#To execute the benchmarks you will need python2.7 and python-tk

> sudo apt install python2.7 python-tk

#Additionally, you need matplotlib in python2.7

> pip install --user matplotlib

#Para revisar que matplotlib esté instalado correctamente, para python2.7
#Abrir python 2.7
#si no lo reconoce como instalado al intentar importarlo, cerrar python
> python2.7
>>import matplotlib
>>exit()
#e instalar matplotlib de la forma:
> sudo apt-get install python-matplotlib

#Para visualizar lar gráficas, luego de hacer "make"
#Situarse en la carpeta de benchmarks dentro de build
> cd /build/benchmarks

#y correr el comando:
> ./benchmark -t RootFinders -r detailed

#Las gráficas aparecerán una a una, cuando se cierra la primera se abre la segunda, y asíi sucesivamente 
#en donde las primeras cuatro gráficas corresponde a las 4 funciones para variables de tipo float
#y en donde las últimas cuatro gráficas corresponden a las 4 funciones para variables de tipo double

