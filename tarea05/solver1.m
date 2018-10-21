function ans = solver1()
  powFactor = (3:10); %Creates an array with the factor of the power 
                     %ej: 3 = 1/(2^3) = 1/8 or 10 = 1/(2^10) = 1/1024
  factors = 1./(2.^powFactor); %Creates the array of the factors
  
  %Create the parameter for the values that are going to be used
  f = @function1;
  xi = 0;
  xf = 1;
  y0 = 1;

  %We calculate the vector solution for each step
  %Numbers start at 3 because that the value of the powFactor
  %Of the step used
  [x3,y3] = RungeKutta_4(f,xi,xf,y0,factors(1));
  [x4,y4] = RungeKutta_4(f,xi,xf,y0,factors(2));
  [x5,y5] = RungeKutta_4(f,xi,xf,y0,factors(3));
  [x6,y6] = RungeKutta_4(f,xi,xf,y0,factors(4));
  [x7,y7] = RungeKutta_4(f,xi,xf,y0,factors(5));
  [x8,y8] = RungeKutta_4(f,xi,xf,y0,factors(6));
  [x9,y9] = RungeKutta_4(f,xi,xf,y0,factors(7));
  [x10,y10] = RungeKutta_4(f,xi,xf,y0,factors(8));
  
  %Graphicate the result of each vector
  figure(1);
  plot(x3,y3,"k;1/8 Factor;",
  x4,y4,"r;1/16 Factor;",
  x5,y5,"g;1/32 Factor;",
  x6,y6,"b;1/64 Factor;",
  x7,y7,"y;1/128 Factor;",
  x8,y8,"m;1/256 Factor;",
  x9,y9,"c;1/512 Factor;",
  x10,y10,"k;1/1024 Factor;");
  
  %Obatin the error of the last point of each vector solution
  err = zeros(1,8);
  err(1) = abs(2-y3(length(y3)));
  err(2) = abs(2-y4(length(y4)));
  err(3) = abs(2-y5(length(y5)));
  err(4) = abs(2-y6(length(y6)));
  err(5) = abs(2-y7(length(y7)));
  err(6) = abs(2-y8(length(y8)));
  err(7) = abs(2-y9(length(y9)));
  err(8) = abs(2-y10(length(y10)));
  
  %Graphicate the error
  figure(2);
  semilogy(factors,err,"-+k");
endfunction
