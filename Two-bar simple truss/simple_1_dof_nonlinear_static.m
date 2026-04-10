
clear all
close all

% opsætning af parametre

a = 1; %[m]
h = 1; %[m]

L0 = sqrt(a^2 + h^2); %[m]

EA1 = 1; %[N] material property of the left truss
alpha = 1; % scaling factor 

e0 = (a-L0)/L0;
N0 = -EA1*e0;

% calculating of points based on the analytical solution

n = 100; % number of iterations

u = zeros(n+1,1); % displacement vector
P = zeros(n+1,1); % load vector

xmax = 2.5; % maximum relative displacement looked at (xmax = u_max/h)

for i = 1:n+1
  u1 = xmax*h*(i-1)/n; % displacement increments
  L = sqrt(a^2 + (h-u1)^2); % length based on displacement
  e1 = (L-L0)/L0;
  N = EA1*e1;
  P1 = -(1+alpha)*((h-u1)/L)*N;
  u(i) = u1/h;
  P(i) = P1/N0;
end

plot(u,P,'-k','linewidth',1)
title('gemetric non-linear realationship for simple truss structure')
xlabel('relative displacement, u/h')
ylabel('relative load, P/N0')

