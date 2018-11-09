/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @Author: 
 * @Date  : 24.02.2018
 */

#include <boost/program_options.hpp>
#include <boost/type_traits/is_complex.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string_regex.hpp>

#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <cstdlib>
#include <iostream>

#include <AnpiConfig.hpp>
#include <Exception.hpp>

namespace po = boost::program_options;

int main(int argc, const char *argv[]) {

	try {
    
    po::options_description desc("Opciones");
	  
	  /*
	  // Parametros a recibir de la interfaz con la consola
		std::vector<int> top, bottom, left, right;
		std::vector<std::string> isolate, file;
		int horizontal, vertical, grid;
		
		//Flags para el aislamiento
		bool explicitTop = false,
		explicitBottom = false,
		explicitLeft = false,
		explicitRight = false;

		// Declare the supported options.
		po::options_description desc("Opciones");
		desc.add_options()
		    ("help", "produce help message")
		    ("top,t", po::value< std::vector<int> >(&top)->multitoken(), "Indica temperatura en borde superior")
		    ("bottom,b", po::value< std::vector<int> >(&bottom)->multitoken(), "Indica temperatura en borde inferior")
		    ("left,l", po::value< std::vector<int> >(&left)->multitoken(), "Indica temperatura en borde izquierdo")
		    ("right,r", po::value< std::vector<int> >(&right)->multitoken(), "Indica temperatura en borde derecho")
		    ("isolate,i", po::value< std::vector<std::string> >(&isolate)->multitoken(), "Aisla los bordes indicados (t=arriba, b=abajo, l=izquierda, r=derecha)")
		    ("profile", po::value< std::vector<std::string> >(&file), "Nombre de archivo con perfil termico")
		    ("horizontal,h", po::value<int>(&horizontal)->default_value(100), "Numero de pixeles horizontales en la solucion")
		    ("vertical,v", po::value<int>(&vertical)->default_value(100), "Numero de pixeles verticales en la solucion")
		    ("quit-visuals,q", "Desactiva toda forma de visualizacion en caso de estar presente")
		    ("flow,f", "Activa el calculo del flujo de calor")
		    ("grid,g", po::value<int>(&grid)->default_value(-1), "Tamaño de rejilla de visualizacion para flujo de calor")
		;

		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);    

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
      if (isolate.empty()) {
        throw anpi::Exception("Aislamiento usado pero no especificado");
      } else {
        std::cout << "Lado aislado es: ";
        for (auto& element : isolate) {
          std::cout << element << " ";
          switch (*element.c_str()) {
            case 't': {
              explicitTop = true;
              break;
            }
            case 'b': {
              explicitBottom = true;
              break;
            }
            case 'l': {
              explicitLeft = true;
              break;
            }
            case 'r': {
              explicitRight = true;
              break;
            }
            default: {
              throw anpi::Exception("Lado incorrecto para aislar");
            };
          }
        }
        std::cout << std::endl;
      }
		}*/
	} catch (const anpi::Exception &exc) {
		std::cerr << exc.what() << '\n';
	}
}

