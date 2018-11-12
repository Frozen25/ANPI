/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @Author: Andres Ramirez
 * @Date  : 08.11.2018
 */

namespace po = boost::program_options;

int top, bottom, left, right;

// Declare the supported options.
po::options_description desc("Opciones");
desc.add_options()
    ("help", "produce help message")
    ("top,t", po::value<int>(&top)->default_value(-1), "Indica temperatura en borde superior")
    ("bottom,b", po::value<int>(&bottom)->default_value(-1), "Indica temperatura en borde inferior")
    ("left,l", po::value<int>(&left)->default_value(-1), "Indica temperatura en borde izquierdo")
    ("right,r", po::value<int>(&right)->default_value(-1), "Indica temperatura en borde derecho")
    ("isolate,i", po::value< vector<string> >(), "Aisla los bordes indicados (t=arriba, b=abajo, l=izquierda, r=derecha)")
    ("profile", po::value< vector<string> >(), "Nombre de archivo con perfil termico")
    ("horizontal,h", po::value<int>()->default_value(100), "Numero de pixeles horizontales en la solucion")
    ("vertical,v", po::value<int>()->default_value(100), "Numero de pixeles verticales en la solucion")
    ("quit-visuals,q", "Desactiva toda forma de visualizacion en caso de estar presente")
    ("flow,f", "Activa el calculo del flujo de calor")
    ("grid,g", po::value<int>()->default_value(-1), "Tamaño de rejilla de visualizacion para flujo de calor")
;

po::variables_map vm;
po::store(po::parse_command_line(ac, av, desc), vm);
po::notify(vm);    

if (vm.count("help")) {
    cout << desc << "\n";
    return 1;
}

if (vm.count("top")) {
    cout << "Temperatura del borde superior es " 
 << vm["top"].as<int>() << ".\n";
}
if (vm.count("bottom")) {
    cout << "Temperatura del borde inferior es " 
 << vm["bottom"].as<int>() << ".\n";
}
if (vm.count("left")) {
    cout << "Temperatura del borde izquierdo es " 
 << vm["left"].as<int>() << ".\n";
}
if (vm.count("right")) {
    cout << "Temperatura del borde derecho es " 
 << vm["right"].as<int>() << ".\n";
}