% Crisptofer Fernández
% 19/10/18
% RungeKutta_4 de Tarea05

% E: function to evaluate (f), initial point (xi), ending point (xf), 
% value in 'y' when the function is evaluated in 0 (y0), size of the step (h)
% S: vector of the values in x and y
function [x,y] = RungeKutta_4(f,xi,xf,y0,h)
  x = xi:h:xf;
  y = zeros(1,length(x));
  y(1) = y0;
  for n = 2:length(x)
    k1 = f(x(n-1),y(n-1));
    k2 = f(x(n-1)+h/2,y(n-1)+h*k1/2);
    k3 = f(x(n-1)+h/2,y(n-1)+h*k2/2);
    k4 = f(x(n-1)+h,y(n-1)+h*k3);
    y(n) = y(n-1)+(h/6)*(k1+2*k2+2*k3+k4);
  endfor
end