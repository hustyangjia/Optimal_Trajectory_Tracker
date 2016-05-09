function Y = test_ode(t,x)
  if t==0
    global Torque
    Torque = [0; 0]
  end
  global Torque

  X=x;
  Y = zeros(size(X));

  R.m_t = 10;
  R.m_w = 2;
  R.I_ybb = 1;
  R.I_t = 5;
  R.r = 0.05;
  R.d = 0;
  R.b = 0.4;

  qd = [0 0 0]';
  vd = 0.5;
  if t<7
  	wd = 1;
  else
  	wd = -1;
  end

  x = X(1);
  y = X(2);
  theta = X(3);
  v = X(4);
  w = X(5);

  W = reshape(X(6:end),15,1);

  ex = x - qd(1);
  exdot = v*cos(theta);
  ey = y - qd(2);
  eydot = v*sin(theta);
  etheta = theta - qd(3);
  ethetadot = w;

  M_q = diag([R.m_t+2*R.I_ybb/R.r^2, R.m_t*R.d^2+R.I_t+2*R.I_ybb*R.b^2/R.r^2-4*R.m_w*R.d^2]);
  B_q = [1 1; R.b -R.b]/R.r;
  G_q = 0;
  F_q = 0;
  V_m = R.d*w*(R.m_t - 2*R.m_w) * [0 -1; 1 0];
  Velocity = [v; w];
  Velocity_temp = inv(M_q)*(B_q*Torque - V_m*Velocity - G_q - F_q);
  
  ev = v - vd;
  evdot = Velocity_temp(1);
  ew = w - wd;
  ewdot = Velocity_temp(2);

  PhiE = [	ex^2
  			ex*ey
  			ex*etheta
  			ex*ev
  			ex*ew
  			ey^2
  			ey*etheta
  			ey*ev
  			ey*ew
  			etheta^2
  			etheta*ev
  			etheta*ew
  			ev^2
  			ev*ew
  			ew^2 	];

  PhiEd = [	2*ex*exdot 0 0 0 0
      			ey*exdot ex*eydot 0 0 0 
      			etheta*exdot 0 ex*ethetadot 0 0
      			ev*exdot 0 0 ex*evdot 0
      			ew*exdot 0 0 0 ex*ewdot
      			0 2*ey*eydot 0 0 0
      			0 etheta*eydot ey*ethetadot 0 0
      			0 ev*eydot 0 ey*evdot 0
      			0 ew*eydot 0 0 ey*ewdot
      			0 0 2*etheta*ethetadot 0 0
      			0 0 ev*ethetadot etheta*evdot 0
      			0 0 ew*ethetadot 0 etheta*ewdot
      			0 0 0 2*ev*evdot 0
      			0 0 0 ew*evdot ev*ewdot
      			0 0 0 0 2*ew*ewdot             ];

  alpha_1 = 1;
  alpha_2 = 0.1;

  E = [ex ey etheta ev ew]';
  Q = E'*E;
  R = eye(4);
  f_e = zeros(5,1);
  GX = zeros(5,4);
  GX(1:3,1:2) = [cos(theta) 0; sin(theta) 0; 0 1];
  GX(4:5,3:4) = inv(M_q)*B_q;
  GG = GX*inv(R)*GX';

  disp(R)
  disp(GX)
  disp(PhiE)
  disp(PhiEd)
  disp(W)

  u = - 1/2 * inv(R) * GX' * PhiEd' * W;
  Torque = u(3:4);

  V_m = 0;
  Velocity_temp = inv(M_q)*(B_q*Torque - V_m*Velocity - G_q - F_q);
  v_d = Velocity_temp(1);
  w_d = Velocity_temp(2);
  
  x_d = v*cos(theta);
  y_d = v*sin(theta);
  theta_d = w;

  Sigma_hat = PhiEd*(f_e + GX*u);
  W1 = -alpha_1*Sigma_hat * (Q + W'*PhiEd*(f_e + GG*PhiEd'*W/4))/(Sigma_hat'*Sigma_hat+1)^2;
  Wd = W1;
  
  Y(1) = x_d;
  Y(2) = y_d;
  Y(3) = theta_d;
  Y(4) = v_d;
  Y(5) = w_d;
  Y(6:end) = Wd(:);
end