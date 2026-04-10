%function Non_linear_CDM
clear all
close all

% Non-linear 1-dof system

%% initiation of known parameters

a_L = 1; %[m] width of bar
h = 1; %[m] height of bar

L0 = sqrt(a_L^2 + h^2); %[m] undeformed length of bar

E = 210e9; %[Pa]
A = (pi/4)*(50e-3)^2; %[m^2]
EA = E*A; %[N] material property of the truss'

k_lin = 2*(EA/L0)*(h/L0)^2;
zeta = 0;%0.04;
m = 30; %[kg]
c_lin = 2*sqrt(k_lin*m)*zeta;
omega_n = sqrt(k_lin/m);
omega_d = omega_n*sqrt(1-zeta^2);
R0 = 0;%10e6; %[N]
r0 = R0/m;
omega_e = 1.1*omega_n;
chi = atan2((2*zeta*omega_n*omega_e),(omega_n^2 - omega_e^2));

% switches 
saver = 1;

%% Initial conditions

u0 = 0.3*h;
du0 = 0;
ddu0_lin = r0 - 2*zeta*omega_n*du0 - omega_n^2*u0;
ddu0_non = (1/m)*(R0 - c_lin*du0 - Q(E,A,a_L,h,L0,u0));
ddu0_non_Lagr = (1/m)*(R0 - c_lin*du0 - Q_L(E,A,h,L0,u0));

% disp('error inducing initial values')
% u01_error = 1;
% u02_error = 2;
% u03_error = 3;


%% Linear solution coefficients

H = r0/(sqrt((omega_n^2 - omega_e^2)^2 + (2*zeta*omega_n*omega_e)^2));
A_lin = u0 - H*cos(-chi);
B_lin = (du0 + zeta*omega_n*A_lin + H*omega_e*sin(-chi))/omega_d;


%% Set-up of iteration parameters
T = 2*pi/omega_n; %Period with respect to the undamped linear solution
Td = T*(1/sqrt(1-zeta^2));
%dt = T/10e5; %Time integration step
%dt_d = Td/10e5; %Time integration step

%N = round(20*T/dt); % number of iterations
if zeta > 0
  dt = Td/10e5;
  t_CDM = 0:dt:20*Td;
  N = length(t_CDM);
else
  dt = T/10e5;
  t_CDM = 0:dt:20*T;
  N = length(t_CDM);
end
Dt = zeros(N+1,1); % Vector of time points

%% Non-linear Central difference method (CDM)
u_non = zeros(N+2,1); % Non-linear displacement
u_non_Lagr = zeros(N+2,1);
du_non = zeros(N+1,1);
ddu_non = zeros(N+1,1);
R = zeros(N+1,1);
Q_ut = zeros(N+1,1);
Q_ut_Lagr = zeros(N+1,1);
dE_mec = zeros(N+1,1);

a_non = (1/dt^2)*m - (1/(2*dt))*c_lin;
b_non = (2/dt^2)*m;
Z_non = ((1/dt^2)*m + (1/(2*dt))*c_lin);

u_before = u0 - du0*dt + (1/2)*ddu0_non*dt^2;
u_non(1) = u_before;
u_non(2) = u0;

du_non(1) = du0;
ddu_non(1) = ddu0_non;

u_before_Lagr = u0 - du0*dt + (1/2)*ddu0_non*dt^2;
u_non_Lagr(1) = u_before_Lagr;
u_non_Lagr(2) = u0;

R(1) = R0;
Q_ut(1) = Q(E,A,a_L,h,L0,u_non(2));
dE_mec(1) = dEdt(m,E,A,a_L,h,L0,ddu_non(1),du_non(1),u_non(2)); 

Q_ut_Lagr(1) = Q_L(E,A,h,L0,u_non_Lagr(2));

for i = 2:N+1
  t = (i-1)*dt;
  Dt(i) = t;

  R(i) = R0*cos(omega_e*t); 

  u_non(i+1) = (Z_non^-1)*(R(i-1) - Q(E,A,a_L,h,L0,u_non(i)) + b_non*u_non(i) - a_non*u_non(i-1));

  du_non(i) = (u_non(i+1) - u_non(i-1))/(2*dt);
  ddu_non(i) = (u_non(i+1) - 2*u_non(i) + u_non(i-1))/dt^2;

  Q_ut(i) = Q(E,A,a_L,h,L0,u_non(i+1));
  dE_mec(i) = dEdt(m,E,A,a_L,h,L0,ddu_non(i),du_non(i),u_non(i+1));

  u_non_Lagr(i+1) = (Z_non^-1)*(R(i-1) - Q_L(E,A,h,L0,u_non_Lagr(i)) + b_non*u_non_Lagr(i) - a_non*u_non_Lagr(i-1));

  Q_ut_Lagr(i) = Q_L(E,A,h,L0,u_non_Lagr(i+1));

end 

%% Linear central difference method

u_lin = zeros(N+2,1);
du_lin = zeros(N+1,1);
ddu_lin = zeros(N+1,1);

u_lin_before = u0 - du0*dt + (1/2)*ddu0_lin*dt^2;
u_lin(1) = u_lin_before;
u_lin(2) = u0;
du_lin(1) = du0;
ddu_lin(1) = ddu0_lin;

a_lin = (1/dt^2)*m - (1/(2*dt))*c_lin;
b_lin = (2/dt^2)*m - k_lin;
Z_lin = ((1/dt^2)*m + (1/(2*dt))*c_lin);


for i = 2:N+1
  u_lin(i+1) = (Z_lin^-1)*(R(i) + b_lin*u_lin(i) - a_lin*u_lin(i-1));
  du_lin(i) = (u_lin(i+1) - u_lin(i-1))/(2*dt);
  ddu_lin(i) = (u_lin(i+1) - 2*u_lin(i) + u_lin(i-1))/dt^2;
end

%% Linear solution

u_lin_ana = zeros(N+1,1); % Linear displacement
u_lin_ana(1) = u0;
u_ana_beats = zeros(N+1,1);
u_ana_beats(1) = u0;

for i = 2:N+1
  t = (i-1)*dt;
  u_lin_ana(i) = exp(-zeta*omega_n*t)*(A_lin*cos(omega_d*t)+B_lin*sin(omega_d*t)) + H*cos(omega_e*t-chi);
end

if R0 > 0 && zeta == 0
  for i = 2:N
    t = (i-1)*dt;
    u_ana_beats(i) = (2*r0/(omega_n^2 - omega_e^2))*sin(((omega_n - omega_e)/2)*t)*sin(((omega_n + omega_e)/2)*t);
  end
end

u_ana_beats_plot = u_ana_beats/max(u_lin_ana);
  
%% Energy tracking

% Non-linear Mechanical energy of the system
EnergyQ_non = zeros(N,1);
EnergyQ_non = Eq(E,A,a_L,h,L0,u_non(2:end));

EnergyK_non = zeros(N,1);
EnergyK_non = Ek(m,du_non);

E_mec_non = EnergyK_non + EnergyQ_non;

dE_mec_CR_non = zeros(N,1);
dE_mec_CR_non = dEcr(c_lin,R,du_non);

% Linear mechanical energy of the system
EnergyQ_lin = zeros(N,1);
EnergyQ_lin = Eql(k_lin,u_lin(2:end));

EnergyK_lin = zeros(N,1);
EnergyK_lin = Ek(m,du_lin);

E_mec_lin = EnergyK_lin + EnergyQ_lin;

% Energy consideration
E_el0 = Eq(E,A,a_L,h,L0,0); %- Eq(E,A,a_L,h,L0,0); %Total mechanical energy of the system at t=0 if du0=0
u_non_zeros = find(abs(u_non(2:end)) < 1e-5);
du_non_zeros = find(abs(du_non) < 1e-5);
length(u_non_zeros);
length(du_non_zeros);
E_kin = Ek(m,du_non(u_non_zeros)); %Total mechanical energy of the system for every t where u(t)=0
E_el =  Eq(E,A,a_L,h,L0,u_non(du_non_zeros)); %- Eq(E,A,a_L,h,L0,0);%Total mechanical energy of the system for every t where du(t)=0


%% Plots

Dt_plot = Dt/T;


% Excitation force plot

% plot of the excitation force
% figure()
% plot(Dt_plot,R/max(R))
% xlabel('Relative time, t/T [s]')
% ylabel('Relative load, R/R_{max} [N]')

% Plots relating to the displacements

u_non_plot = u_non(2:end)/abs(max(u_non(2:end)));
u_non_Lagr_plot = u_non_Lagr(2:end)/abs(max(u_non(2:end)));
u_lin_plot = u_lin(2:end)/abs(max(u_non(2:end)));
u_lin_ana_plot = u_lin_ana/abs(max(u_non(2:end)));


figure()
plot(Dt_plot,u_non_plot,':k.')
hold on
plot(Dt_plot,u_lin_ana_plot,'-k')
hold on
plot(Dt_plot,u_lin_plot,'--k')
xlabel('Relative time, t/T [s]')
ylabel('Relative displacement, u/u_{lin,max} [m]')
title('Displacement over time')
legend({'CDM non-linear, u_{non}','linear Analytisk, u_{lin}','CDM linear'},'Location','northeast')

hold off
figure()
plot(Dt_plot,u_non_plot,':k.')
hold on
plot(Dt_plot,u_non_Lagr_plot,'-k')
xlabel('Relative time, t/T')
ylabel('Relative displacement, u/u_{max}')
title('Displacement over time')
legend({'CDM Non-linear, u_{non}','CDM Non-linear Lagrange, u_{non,lagr}'},'Location','northeast')

hold off
figure()
plot(Dt_plot,u_non_Lagr_plot,':k.')
hold on
plot(Dt_plot,u_lin_ana_plot,'-k')
hold on
plot(Dt_plot,u_lin_plot,'--k')
xlabel('Relative time, t/T [s]')
ylabel('Relative displacement, u/u_{lin,max} [m]')
title('Displacement over time')
legend({'CDM non-linear Lagrange, u_{non}','linear Analytisk, u_{lin}','CDM linear'},'Location','northeast')

% plots relating to the non-linear internal restoring force of the system

hold off
figure()
plot(u_non(2:end),Q_ut,':k.')
xlabel('displacement, u_{non}')
ylabel('internal restoring force')
title('internal force in dynamic motion')

figure()
plot(u_non(2:end),Q_ut,'-k')
hold on
plot(u_non_Lagr(2:end),Q_ut_Lagr,'--k');
xlabel('displacement')
ylabel('Internal retoring force')
title('Internal force in dynamic motion with lagrange')
legend({'Analytical force, Q_{ut}','Lagrange force, Q_{ut,L}'},'Location','northeast')

hold off

% plots relating to the energy of the system

E_mec_non_plot = E_mec_non/max(E_mec_non);
E_mec_lin_plot = E_mec_lin/max(E_mec_lin);

figure()
plot(Dt_plot,E_mec_non_plot,':k.')
axis equal
xlabel('Relative time, t/T')
ylabel('Relative energy, E/E_{max}')
title('Mechanical energy of the non-linear system over time')

figure()
plot(Dt_plot,E_mec_lin_plot,':k.')
axis equal 
xlabel('Relative time, t/T')
ylabel('Relative energy, E/E_{max}')
title('Mechanical energy of the linear system over time')

figure()
plot(Dt_plot,dE_mec,':k.')
xlabel('time, t [s]')
ylabel('Energy/time, dE/dt [J/s]')
title('Change in mechanical energy over time')

figure()
% hold on
plot(Dt_plot,dE_mec_CR_non,'-k')
xlabel('time, t [s]')
ylabel('Energy/time, dE/dt [J/s]')
title('Change in mechanical energy over time (force and damping term)')

%% Test of the stiffness function

u_test = linspace(-1,2.5,1000);
Q_u = zeros(length(u_test),1);
Q_u_Lagr = zeros(length(u_test),1);
Q_u_lin = zeros(length(u_test),1);
for i = 1:length(u_test)
  Q_u_lin(i) = k_lin*u_test(i);
  Q_u(i) = Q(E,A,a_L,h,L0,u_test(i));
  Q_u_Lagr(i) = Q_L(E,A,h,L0,u_test(i));
end
u_test = u_test'/h;

figure()
plot(u_test,Q_u,':k.')
hold on
plot(u_test,Q_u_lin,'--k')
hold on
plot(u_test,Q_u_Lagr,'-k')
xlabel('u/h [-]')
ylabel('force [N]')
title('Non-linear internal restoring force vs linear force')
legend({'Non-linear analytical Q, Q_{non}','linear Analytisk Q, Q_{lin}','Non-linear Lagrange Q, Q_L'},'Location','northeast')

%% Test of the internal non-linear energy function

u_test1 = linspace(-1.5,2.5,300);
EnergyQ_test = zeros(length(u_test1),1);
EnergyQL_test = zeros(length(u_test1),1);
EnergyQlin_test = zeros(length(u_test1),1);
for i = 1:length(u_test1)
  EnergyQ_test(i) = Eq(E,A,a_L,h,L0,u_test1(i)); %- Eq(E,A,a_L,h,L0,0);
  EnergyQL_test(i) = Eq_L(E,A,h,L0,u_test1(i));
  EnergyQlin_test(i) = Eql(k_lin,u_test1(i));
end
u_test1 = u_test1';

figure()
plot(u_test1,EnergyQ_test,':k.')
hold on
plot(u_test1,EnergyQL_test,'-k')
hold on
plot(u_test1,EnergyQlin_test,'--k')
xlabel('u [m]')
ylabel('Energy [J]')
title('Non-linear internal potential energy')

%% The eigenvalue problem

Omega = zeros(N,1);
Omega_ee = zeros(N+1,1);
Omega_n = omega_n*ones(N,1);
Omega_e = omega_e*ones(N,1);

for i = 1:N
  Omega(i) = kt(E,A,a_L,h,L0,u_non(i))/m;
end

omega_non = sqrt(Omega);
Omega_mean = mean(Omega);
omega_non_mean = mean(omega_non);

Omega_non_mean = omega_non_mean*ones(N,1);

u_ana_beats_non = zeros(N,1);
u_ana_beats_non(1) = u0;

if R0 > 0 && zeta == 0
  for i = 2:N
    t = i*dt;
    u_ana_beats_non(i) = (2*r0/(Omega_mean - omega_e^2))*sin(((omega_non_mean - omega_e)/2)*t)*sin(((omega_non_mean + omega_e)/2)*t);
  end
end

u_ana_beats_non_plot = u_ana_beats_non/max(u_lin_ana);

% hold off
% 
% figure()
% plot(Dt_plot,omega_non,'--k')
% hold on
% plot(Dt_plot,Omega_n,':k.')
% hold on
% plot(Dt_plot,Omega_e,'-k')
% hold on
% plot(Dt_plot,Omega_non_mean,':k.')
% xlabel('time, t [s]')
% ylabel('frequencies, omega [Hz]')
% title('non-linear frequencies vs linear frequency')
% 
% hold off
% 
% figure()
% plot(Dt_plot,u_ana_beats_plot,'-k')
% hold on
% plot(Dt_plot,u_ana_beats_non_plot,'-k')
% xlabel('Relative time ,t/T')
% ylabel('Relative displacement, u/u_{max}')
% title('non-linear beats vs linear beats')

%% Saving

if saver == 1
  if R0 == 0 && zeta == 0
    save('1dof_CDM_SimpleHarmonic')
  elseif R0 == 0 && zeta > 0
    save('1dof_CDM_DampedHarmonic')
  elseif R0 > 0 && zeta > 0
    save('1dof_CDM_ForcedDampedHarmonic')
  end
end



%% Convinient functions

% Internal non-linear restoring forces (non-linear geometric)
function Q = Q(E,A,a_L,h,L0,u)
  EA = E*A;
  L = sqrt(a_L^2 + (h-u)^2);
  Q = -2*EA*((h-u)/L)*(L/L0-1);
end

% Tangent stiffnes coefficient (non-linear geometric)
function kt = kt(E,A,a_L,h,L0,u)
  EA = E*A;
  L = sqrt(a_L^2 + (h-u)^2);
  NN = (EA/L0)*(L-L0);
  dN = -(EA/L0)*(h-u)/L;
  dhu = (-1/L)*(1-((h-u)/L)^2);
  kt = -2*(dN*((h-u)/L) - NN*dhu);
end

% Internal force with the lagrange strain measure
function Q_L = Q_L(E,A,h,L0,u)
  EA = E*A;
  e = (-(h/L0) + (1/2)*(u/L0))*(u/L0);
  N = EA*e;
  Q_L = -2*N*(h-u)/L0;
end

% Internal non-linear energy (non-linear geomeric)
% function Eq = Eq(E,A,a_L,h,L0,u)
%   EA = E*A;
%   L = sqrt(a_L^2 + (h-u).^2);
%   Eq = -2*EA*((u/L0).*(h-u/2) + L);
% end

function Eq = Eq(E,A,a_L,h,L0,u)
  EA = E*A;
  L = sqrt(a_L^2 + (h-u).^2);
  Eq = 2*EA*(L0 - L + (u-2*h).*u/(2*L0));
end

% Internal nonlinear energy (Lagrange)
function Eq_L = Eq_L(E,A,h,L0,u)
  EA = E*A;
  Eq_L = (EA/L0)*(u.^2/4 - h*u + h^2)*(u/L0)^2;
end

% Internal linear potential energy
function Eql = Eql(k_lin,u)
  Eql = (1/2)*k_lin*u.^2;
end

% kinetic energy
function Ek = Ek(m,du)
  Ek = (1/2)*m*du.^2;
end

% Change in mechanical energy over time
function dEdt = dEdt(m,E,A,a_L,h,L0,ddu,du,u)
  EA = E*A;
  L = sqrt(a_L^2 + (h-u).^2);
  dEdt = m.*ddu.*du + 2*EA.*du.*(h-u).*((1/L) - (1/L0));
end

% change in mechanical energy over time represented by the damping term and
% the excitation force
function dEc = dEcr(c,R,du)
  dEc = R.*du - c.*du.^2;
end

%end























