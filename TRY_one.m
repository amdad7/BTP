% Finite element calculation - Forced Vibration - Cantilever Beam

clear;
clc;

global f2;
n = 4;
%initial conditions
IC=zeros(4*(n),1);

%Time SPan(seconds)
t0=0; tf=2;
tspan=[t0 ,tf];

%Numerical integration
[time,state_values]=ode45(@g,tspan,IC);
x=state_values(:,2);
vel=state_values(:,10);
theta=state_values(:,6);
ang=state_values(:,14);

%FFT
Y = fft(x);

%Plot results
figure(1), clf
plot(time,x), xlabel('time (s)'),ylabel('displacement (m)')
title('X Displacement vs. Time')

figure(2), clf
plot(time, theta), xlabel('time(s)'), ylabel('angular displacement (rad)')
title('\theta angular dispalacement vs. time')

figure(3), clf
plot(time, vel), xlabel('time(s)'), ylabel('velocity (m/s)')
title('V Velocity vs. Time')

figure(4), clf
plot(time, ang), xlabel('time(s)'), ylabel('angular velocity (m/s)')
title('\omega Angular Velocity vs. Time')

figure(5)
plot(Y)
f2

