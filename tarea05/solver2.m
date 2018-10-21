function [ans45,ans23,ans4] = solver2()
  %Create the parameter for the values that are going to be used
  f = @function2;
  trange = [0,200];
  y0 = 5;
  
  %Initialize the return value in the form (steps,time)
  ans45 = zeros(1,2);
  ans23 = zeros(1,2);
  ans4 = zeros(1,2);
  
  %The number corresponds to the ode method used and measures the time it took
  %To execute each method
  tic();
  [x45,y45] = ode45(f,trange,y0);
  ans45(2) = toc();
  
  tic();
  [x23,y23] = ode23(f,trange,y0);
  ans23(2) = toc();
  
  tic();
  [x4,y4] = RungeKutta_4(f,trange(1),trange(2),y0,0.2);
  ans4(2) = toc();
  %By using 0.2 we are going to use at least 1000 steps
  
  %Returns the amount of steps used by each method
  ans45(1) = length(x45);
  ans23(1) = length(x23);
  ans4(1) = length(x4);
  
  %Graphicate the result of each vector
  figure(1);
  plot(x45(40:length(x45)),y45(40:length(y45)),"b;ode45;",
  x23(54:length(x23)),y23(54:length(y23)),"r;ode23;",
  x4(501:length(x4)),y4(501:length(y4)),"g;Runge Kutta 4;");
endfunction
