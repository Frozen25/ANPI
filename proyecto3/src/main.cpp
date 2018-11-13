/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @Author: Andres Ramirez-Quiros
 * @Date  : 07.11.2018
 */

#include <boost/program_options.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>

#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <cstdlib>
#include <iostream>

#include <AnpiConfig.hpp>
#include <Exception.hpp>
#include <fstream>

#include <ThermalPlate.hpp>

namespace po = boost::program_options;

void fileParser(std::string &file, std::vector<float> &top, std::vector<float> &bottom,
                std::vector<float> &left, std::vector<float> &right) {
  std::ifstream inFile(file.c_str());
  
  if (!inFile) {
    throw anpi::Exception("No se pudo abrir el archivo");
  } else {
    std::cout << "Archivo contiene la siguiente informacion: " << std::endl;
  }
  
  std::vector<float> *vectorToUse; // Puntero a vector en el cual escribir. Depende del lado asignado
  
  std::string line;
  getline(inFile, line);
  
  do {
    std::vector<std::string> temperature;
    boost::split(temperature, line, [](char c){return c == ' ';});
  
    switch ((temperature[0].c_str())[0]) {
      case 't': {
        vectorToUse = &top;
        std::cout << "Top: ";
        break;
      }
      case 'b': {
        vectorToUse = &bottom;
        std::cout << "Bottom: ";
        break;
      }
      case 'l': {
        vectorToUse = &left;
        std::cout << "Left: ";
        break;
      }
      case 'r': {
        vectorToUse = &right;
        std::cout << "Right: ";
        break;
      }
      default: {
        throw anpi::Exception("Archivo vacio o caracter invalido");
      }
    }
  
    if (vectorToUse->empty()) {
      for (size_t i = 2; i < temperature.size(); ++i) { // Empieza en 2 pues ignora la palabra de lado y el signo =
        try {
          vectorToUse->push_back(std::stoi(temperature[i]));
          std::cout << temperature[i] << " ";
        } catch (const std::invalid_argument& ia) {
            std::cerr << "Argumento invalido: " << ia.what() << '\n';
          }
      }
    } else {
      std::cout << "temperatura no será utilizada";
    }
    std::cout << std::endl;
    
    getline(inFile, line);
  } while (!line.empty());
  
  inFile.close();
}

int main(int argc, const char *argv[]) {

	try {
	  // Parametros a recibir de la interfaz con la consola
		std::vector<float> top, bottom, left, right;  //En caso que no se asignen temperaturas, se utilizará empty() para ver si es aislado
		std::vector<std::string> isolate;
		std::string file;
		int horizontal, vertical, grid;
		
		// Flags para el aislamiento
		// Los ExplicitFlags es para cuando se asigna espeficiamente que un lado esta aislado
		// Esto para que la lectura del archivo no sobreescriba lo del comando
		bool explicitTop = false,
		explicitBottom = false,
		explicitLeft = false,
		explicitRight = false,
		quitVisuals = false,
		activateFlow = false;

		// Declara las opciones soportadas
		po::options_description desc("Opciones");
		desc.add_options()
		    ("help", "Imprime esta lista de opciones")
		    ("top,t", po::value< std::vector<float> >(&top)->multitoken(), "Indica temperatura en borde superior")
		    ("bottom,b", po::value< std::vector<float> >(&bottom)->multitoken(), "Indica temperatura en borde inferior")
		    ("left,l", po::value< std::vector<float> >(&left)->multitoken(), "Indica temperatura en borde izquierdo")
		    ("right,r", po::value< std::vector<float> >(&right)->multitoken(), "Indica temperatura en borde derecho")
		    ("isolate,i", po::value< std::vector<std::string> >(&isolate)->multitoken(), "Aisla los bordes indicados (t=arriba, b=abajo, l=izquierda, r=derecha)")
		    ("profile,p", po::value< std::string >(&file), "Nombre de archivo con perfil termico")
		    ("horizontal,h", po::value<int>(&horizontal)->default_value(100), "Numero de pixeles horizontales en la solucion")
		    ("vertical,v", po::value<int>(&vertical)->default_value(100), "Numero de pixeles verticales en la solucion")
		    ("quit-visuals,q", "Desactiva toda forma de visualizacion en caso de estar presente")
		    ("flow,f", "Activa el calculo del flujo de calor")
		    ("grid,g", po::value<int>(&grid)->default_value(5), "Tamaño de rejilla de visualizacion para flujo de calor")
		;

		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);
    
    //Creacion de archivo para guardar las flags
    std::ofstream flags;
    flags.open ("flags.txt");

		if (vm.count("help")) {
		  std::cout << desc << "\n";
		  return 1;
		}
		
		if (vm.count("top")) {
		  std::cout << "Temperatura del borde superior es: ";
      for (auto& element : top) {
        std::cout << element << " ";
      }
      std::cout << std::endl;
		}
		
		if (vm.count("bottom")) {
		  std::cout << "Temperatura del borde inferior es: ";
      for (auto& element : bottom) {
        std::cout << element << " ";
      }
      std::cout << std::endl;
		}
		
		if (vm.count("left")) {
		  std::cout << "Temperatura del borde izquierdo es: ";
      for (auto& element : left) {
        std::cout << element << " ";
      }
      std::cout << std::endl;
		}
		
		if (vm.count("right")) {
		  std::cout << "Temperatura del borde derecho es: ";
      for (auto& element : right) {
        std::cout << element << " ";
      }
      std::cout << std::endl;
		}
		
		if (vm.count("isolate")) {
      std::cout << "Lado aislado es: ";
      for (auto& element : isolate) {
        std::cout << element << " ";
        switch (*element.c_str()) {
          case 't': { // Lado superior (top) aislado
            explicitTop = true;
            break;
          }
          case 'b': { // Lado inferior (bottom) aislado
            explicitBottom = true;
            break;
          }
          case 'l': { // Lado izquierdo (left) aislado
            explicitLeft = true;
            break;
          }
          case 'r': { // Lado derecho (right) aislado
            explicitRight = true;
            break;
          }
          default: {
            std::cout << std::endl;
            throw anpi::Exception("Lado incorrecto para aislar");
          };
        }
      }
      std::cout << std::endl;
		}
    
    if (vm.count("profile")) {
      std::cout << "Usando archivo de perfil" << std::endl;
      fileParser(file, top, bottom, left, right);
    }
		
    if (vm.count("quit-visuals")) {
      quitVisuals = true;
      std::cout << "Opcion sin visuales: " << quitVisuals << std::endl;
    }
    
    if (vm.count("flow")) {
      activateFlow = true;
      std::cout << "Opcion para ver el flujo de calor: " << activateFlow << std::endl;
    }
    
    //Flags necesarias para la realizacion del gráfico que se guardan en un archivo aparte
    flags << horizontal << '\n';
    flags << vertical << '\n';
    flags << quitVisuals << '\n';
    flags << activateFlow << '\n';
    flags << grid << '\n';
    
    flags.close();
    
    anpi::ThermalPlate thermalPlate((explicitTop||top.empty()),(explicitBottom||bottom.empty()),
                                    (explicitLeft||left.empty()),(explicitRight||right.empty()),
                                    top, bottom, left, right, activateFlow, quitVisuals,
                                    grid, horizontal, vertical);
    
    int iterations; //Cantidad de iteraciones a realizar según el tamaño deseado de la gráfica
    
    if((horizontal || vertical)>4098){
      iterations = 13;
    } else if((horizontal>2050)||(vertical>2050)){
      iterations = 12;
    } else if((horizontal>1026)||(vertical>1026)){
      iterations = 11;
    } else if((horizontal>514)||(vertical>514)){
      iterations = 10;
    } else if((horizontal>258)||(vertical>258)){
      iterations = 9;
    } else if((horizontal>130)||(vertical>130)){
      iterations = 8;
    } else if((horizontal>66)||(vertical>66)){
      iterations = 7;
    } else if((horizontal>34)||(vertical>34)){
      iterations = 6;
    } else if((horizontal>18)||(vertical>18)){
      iterations = 5;
    } else if((horizontal>10)||(vertical>10)){
      iterations = 4;
    } else if((horizontal>6)||(vertical>6)){
      iterations = 3;
    } else if((horizontal>4)||(vertical>4)){
      iterations = 2;
    } else if((horizontal>3)||(vertical>3)){
      iterations = 1;
    } else {
      iterations = 0;
    }
    thermalPlate.solvePlate( 20.0f , iterations, 1);
    
	} catch (const anpi::Exception &exc) {
		std::cerr << exc.what() << '\n';
	}
}

