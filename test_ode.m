function Y = test_ode(X,T)
  Y = zeros(size(X));

  x = X(1);
  z = X(2);

  u = 1;      % Dummy input

  xd = -x*(pi/2 + atan(5*x)) - (5*x^2)/(2*(1+25*x^2)) + 4x + z;
  zd = 2*x^2 - x + (1 + 0.5*cos(x))*u;

  Y(1) = xd;
  Y(2) = zd;
end
