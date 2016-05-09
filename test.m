clear all
close all
clc

T_end = 10;

W = rand(15,1);

Initial_states = [0.5; -0.5; pi/2; 0; 0; W(:)];
T_span = [0 T_end];

[T,Y] = ode45(@(t,x)  test_ode(t,x), T_span, Initial_states);

plot(Y(:,1),Y(:,2),'x');