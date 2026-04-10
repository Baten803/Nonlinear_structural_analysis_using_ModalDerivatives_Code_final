% Non-linear FEM (Modal derivatives - Lagrange) 

clear all
close all


% Tilføjes damping holdes s-relationerne stadig
% Tilføjes flere mode-shapes vil de samme s-relationer ikke længere holde
% da der er flere termer der nu er med

%% Input data

% Material parameters
E = 210e9; % [Pa] Young's modulus
rho = 15.6; % [kg/m] Density of the element
c_damp = 0; % damping coefficient of the element

% Geometric element parameters
d = 50e-3; % [m] diameter of the cable
A = (pi/4)*(d)^2; % [m^2] Cross sectional area
I = pi*(d)^4/64; % [m^4] second moment of inertia 
a_L = 1; % [m] length of element
L_tot = 20; % [m] Total length of element (use if easier)

% Load parameters
F1 = -1; % [N] Point load
F2 = -3; % [N/m] Line load
e = 0; % [m] Eccentricity
p = rho*9.81;

% Mode parameters
N_phi = 3;
MD_phi = 1;

% Switches
linear_phi = 0;
MD = 0;
V_plot = 1;
Q_plot = 0;
s_relations = 0;


%% Matrices for the element 

% knudekoordinater: X(knude,:)=[x,y] 
Numb = 20; % Number of elements wanted in a straight line

X = zeros(Numb+1,2);

for i = 1:Numb+1
  X(i,:) = [(i-1)*L_tot/Numb 0];
end

nno = size(X,1);    % number of nodes, antal knuder i systemet

% elementer: T(el,:)=[startknude slutknude]

T = zeros(Numb,2);

for i = 1:Numb
  T(i,:) = [i i+1];
end

nel = size(T,1);    % number of elements, antal elementer i systemet

% frihedsgrader, globale dof: D(el,:)=[V1 V2 V3 V4 V5 V6]

D = zeros(nel,6);

for i = 1:nno-1
  D(i,:) = [3*i-2 3*i-1 3*i 3*i+1 3*i+2 3*i+3];
end
 
nd=max(max(D));     % number of dofs, antal frihedsgrader

% materialer: G(el,:)=[E-modul, tværsnitsareal, inertimoment]

G = zeros(nel,5);

for i = 1:nel
  G(i,:) = [E A I rho c_damp];
end

% understøtninger: U(i)=global dof
U(1) = 1;
U(2) = 2;
U(3) = nno*3-2;
U(4) = nno*3-1;

% nodal load (boundary Load): bL(i,:)=[global_dof, magnitude]
bL = []; %empty if no load
%bL = [nno*3-1 -P];
% element load (domain Load): dL(el,:)=[local_direction, magnitude]
 dL_dynamic = zeros(nel,2);
 dL_static = zeros(nel,2);
for i = 11:nel
    dL_dynamic(i,:) = [2 -1];
end 
for i = 1:nel
    dL_static(i,:) = [2 -p];
end

%% plot geometry
plotTop(X,T,nel,nno)

%% Program

% setup displacement vector
V = zeros(nno*3,1);

% setup constant coefficient matrices of the system
M_mass = M(X,T,G,D,nd,nel); % Mass matrix
C_damp = C(X,T,G,D,nd,nel); % Damping matrix
K0 = Kt(X,T,G,nel,nno,V); % standard linear stiffness matrix
R0_dynamic = R(X,T,D,dL_dynamic,bL,nd,nel); % Load amplitude vector
R0_static = R(X,T,D,dL_static,bL,nd,nel);

%% solve linear eigenvalue problem for the EOM
dof = 1:3*nno;                         % index to all dofs 
du = U;                                % index to prescribed dofs
df = setdiff(dof,du);                  % index to free dofs

K0ff = K0(df,df);                          % Free dofs only 
Mff = M_mass(df,df);                          % Free dofs only

[eVdyn,eDdyn] = eig(K0ff,Mff); % solve the eigenvalue problem
omega_sq = sqrt(diag(eDdyn)'); % calculate the natural angular frequencies
[Omega,i] = sort(omega_sq); % Sort the frequencies from lowest to highest
%display(Omega)
eVdyn = eVdyn(:,i); % Sort the modeshapes to the corresponding frequency
                    % i.e from the modeshape belonging to the lowest
                    % frequency to the modeshape belonging to the highest
                    % frequency

phi = zeros(3*nno,size(eVdyn,2));
phi(df,:) = eVdyn; % Here we have added zeros to the non-free dofs such that
                 % we have the full modeshape
Phi = zeros(nd,nd);
Phi(:,df) = phi; % full modal matrix with the zero modal vectors due to boundary conditions
                 % These are normalized with respect to the M-matrix we want it to be the K0-matrix   

m_phi = phi'*M_mass*phi; % making the modal matrices based on the new mode-shapes
c_phi = phi'*C_damp*phi;
k0_phi = phi'*K0*phi;



% Plotting selected modes
Vskala = 1e2;
if linear_phi == 1
  % First mode
  plotDof(X,T,D,phi(:,1),nel,Vskala)            % plot deformed frame

  % Second mode
  plotDof(X,T,D,phi(:,2),nel,Vskala)            % plot deformed frame

  % Third mode
  plotDof(X,T,D,phi(:,3),nel,Vskala)            % plot deformed frame
end


%% solve complete modal derivative eigenvalue problem

% For the length(df) number of mode-shapes in the system
% and H = length(df)*(length(df)-1)/2 number of modal drivatives in the
% system (when V0=0)

H = N_phi*(N_phi+1)/2; % corresponding number of modal derivatives
dPhi = zeros(nd,H);
for i = 1:N_phi
  for j = 1:i
    dKsi = dKs(X,T,G,nel,nno,phi(:,i),V);
    dKsj = dKs(X,T,G,nel,nno,phi(:,j),V);
    alpha0 = phi(df,j)'/m_phi(j,j)*(dKsi(df,df)*phi(df,j) - dKsj(df,df)*phi(df,i));
    dphi = (K0ff - (Omega(i)+Omega(j))^2*Mff)\((alpha0*Mff - dKsi(df,df))*phi(df,j)*2);
    if i==j
      dphi = dphi/2;
    end

    % Map (i,j) to global modal derivative index
    k = i*(i-1)/2 + j;

    dPhi(df,k) = dphi;  
  end
end

if MD == 1
  % % j=k=1
  % dKs11 = dKs(X,T,G,nel,nno,Phi0(:,3),V);
  % alpha10 = phi0(:,1)'/m_phiff(1,1)*(dKs11(df,df)*phi0(:,1) - dKs11(df,df)*phi0(:,1));
  % dphi11 = (K0ff - (Omega(1) + Omega(1))^2*Mff)\(alpha10*Mff - dKs11(df,df))*phi0(:,1)*2;
  % dphi11 = dphi11/2;
  % dPhi11 = zeros(nd,1);
  % dPhi11(df) = dphi11; 
  % 
  dPhi11_plot = dPhi(1:3:end,1)*1e6;
  figure()
  plot(X(:,1),dPhi11_plot,':k.')
  xlabel('x_axis, [m]')
  ylabel('MD 11, [mu*m]')
  title('modal derivative dphi_{11}')
  legend('dphi_{11}','Location','northeast')
  hold off
  % 
  % % j=1, k=2
  % dKs2 = dKs(X,T,G,nel,nno,Phi0(:,4),V);
  % dKs1 = dKs(X,T,G,nel,nno,Phi0(:,3),V);
  % alpha20 = phi0(:,1)'/m_phiff(1,1)*(dKs2(df,df)*phi0(:,1) - dKs1(df,df)*phi0(:,2));
  % dphi21 = (K0ff - (Omega(1) + Omega(2))^2*Mff)\(alpha20*Mff - dKs2(df,df))*phi0(:,1)*2;
  % dPhi21 = zeros(nd,1);
  % dPhi21(df) = dphi21;
  % 
  dPhi21_plot = dPhi(1:3:end,2)*1e6;
  figure()
  plot(X(:,1),dPhi21_plot,':k.')
  xlabel('x_axis, [m]')
  ylabel('MD 21, [mu*m]')
  title('modal derivative dphi_{21}')
  legend('dphi_{21}','Location','northeast')
  hold off
  % 
  % % j=k=2
  % dKs22 = dKs(X,T,G,nel,nno,Phi0(:,4),V);
  % alpha30 = phi0(:,2)'/m_phiff(2,2)*(dKs22(df,df)*phi0(:,2) - dKs22(df,df)*phi0(:,2));
  % dphi22 = (K0ff - (Omega(2) + Omega(2))^2*Mff)\(alpha30*Mff - dKs22(df,df))*phi0(:,2)*2;
  % dphi22 = dphi22/2;
  % dPhi22 = zeros(nd,1);
  % dPhi22(df) = dphi22;
  % 
  dPhi22_plot = dPhi(1:3:end,3)*1e6;
  figure()
  plot(X(:,1),dPhi22_plot,':k.')
  xlabel('x_axis, [m]')
  ylabel('MD 22, [mu*m]')
  title('modal derivative dphi_{22}')
  legend('dphi_{22}','Location','northeast')
  hold off
  % 
  % % j=1, k=3
  % dKs3 = dKs(X,T,G,nel,nno,Phi0(:,5),V);
  % alpha40 = phi0(:,1)'/m_phiff(1,1)*(dKs3(df,df)*phi0(:,1) - dKs1(df,df)*phi0(:,3));
  % dphi31 = (K0ff - (Omega(1) + Omega(3))^2*Mff)\(alpha40*Mff - dKs3(df,df))*phi0(:,1)*2;
  % dPhi31 = zeros(nd,1);
  % dPhi31(df) = dphi31;
  % 
  dPhi31_plot = dPhi(1:3:end,4)*1e6;
  figure()
  plot(X(:,1),dPhi31_plot,':k.')
  xlabel('x_axis, [m]')
  ylabel('MD 31, [mu*m]')
  title('modal derivative dphi_{31}')
  legend('dphi_{31}','Location','northeast')
  hold off
  % 
  % % j=2, k=3
  % alpha50 = phi0(:,2)'/m_phiff(2,2)*(dKs3(df,df)*phi0(:,2) - dKs2(df,df)*phi0(:,3));
  % dphi32 = (K0ff - (Omega(2) + Omega(3))^2*Mff)\(alpha50*Mff - dKs3(df,df))*phi0(:,2)*2;
  % dPhi32 = zeros(nd,1);
  % dPhi32(df) = dphi32;
  % 
  dPhi32_plot = dPhi(1:3:end,5)*1e6;
  figure()
  plot(X(:,1),dPhi32_plot,':k.')
  xlabel('x_axis, [m]')
  ylabel('MD 32, [mu*m]')
  title('modal derivative dphi_{32}')
  legend('dphi_{32}','Location','northeast')
  hold off
  % 
  % % j=k=3
  % alpha60 = phi0(:,3)'/m_phiff(3,3)*(dKs3(df,df)*phi0(:,3) - dKs3(df,df)*phi0(:,3));
  % dphi33 = (K0ff - (Omega(3) + Omega(3))^2*Mff)\(alpha60*Mff - dKs3(df,df))*phi0(:,3)*2;
  % dphi33 = dphi33/2;
  % dPhi33 = zeros(nd,1);
  % dPhi33(df) = dphi33;
  % %dPhi33(df) = -dphi33; % ADD A MINUS SIGN!!!!!!!!!!(not relevant pt)
  % 
  dPhi33_plot = dPhi(1:3:end,6)*1e6;
  figure()
  plot(X(:,1),dPhi33_plot,':k.')
  xlabel('x_axis, [m]')
  ylabel('MD 33, [mu*m]')
  title('modal derivative dphi_{33}')
  legend('dphi_{33}','Location','northeast')
  hold off
  
  % dPhi = [dPhi11 dPhi21 dPhi31 dPhi32 dPhi33]; % putting the MD's into a matrix with all dofs
  % dphi = [dPhi11(df) dPhi21(df) dPhi31(df) dPhi32(df) dPhi33(df)]; % putting the MD's into a matrix with only the free dofs
  % dphi_lincomb = [dPhi11(df) dPhi21(df) dPhi22(df) dPhi31(df) dPhi32(df) dPhi33(df)];
  
  % Test of relation
  
  dphi22_test = 3/8*dPhi(df,4) - 1/4*dPhi(df,1);
  dPhi22_test = zeros(nd,1);
  dPhi22_test(df) = dphi22_test;
  dPhi22_test_plot_x = dPhi22_test(1:3:end)*1e6; % Axiel retning
  dPhi22_test_plot_y = dPhi22_test(2:3:end)*1e10;% Vertical retning
  dPhi22_test_plot_rot = dPhi22_test(3:3:end)*1e10; % rotation
  dphi22_error = dPhi(df,3) - dphi22_test;
  
  
  
  figure()
  plot(X(:,1),dPhi22_plot,'k.')
  hold on
  plot(X(:,1),dPhi22_test_plot_x,'--k.')
  xlabel('x-axis, [m]')
  ylabel('MD 22, dphi22 [m]')
  title('dphi22 vs. 3/8*dphi31-1/4*dphi11')
  legend({'dphi_{22}','3/8*dphi_{31} - 1/4*dphi_{11}'},'Location','northeast')
  hold off
end

%% setting up the modal basis matrix and finding dependent modes

% Initialize the PHI-matrix containing the modes and the modal derivatives

PHI = [phi(:,1:N_phi) dPhi(:,1:end)];
PHi = [phi(df,1:N_phi) dPhi(df,1:end)];

orig_idx = 1:size(PHi,2);
dependent_idx = [];
Beta = zeros(size(PHi,2)-1,size(PHi,2));

if MD_phi == 1
  i = size(phi(:,1:N_phi),2) + 1;
  
  while i <= size(PHI, 2)
      PHI_test = PHI;
      PHI_test(:,i) = [];
  
      beta = (PHI_test' * PHI_test) \ (PHI_test' * PHI(:,i));
      Top = PHI(:,i) - PHI_test * beta;
      e_error = norm(Top) / norm(PHI(:,i)) * 100;
  
      if e_error <= 1
          dependent_idx(end+1) = orig_idx(i);
          PHI(:,i) = [];
          orig_idx(i) = [];
          Beta(:,dependent_idx) = beta; 
      else
          i = i + 1;
      end
  end
  
  PHi = PHI(df,:);
else
  PHI = [phi(:,1:N_phi)];
  PHi = [phi(df,1:N_phi)];
end


%% Set-up of iteration parameters
% T_p = 2*pi/omega_n1; %Period with respect to the undamped linear solution
% dt = T_p/10e4; %Time integration step
% 
% N = round(10*T_p/dt); % number of iterations
% Dt = zeros(N,1); % Vector of time points

dt = 1e-4;
t_CDM = 0:dt:10;
N = length(t_CDM);

%% solve for independent s-functions using CDM

omega_n1 = Omega(1);
omega_e = 1.1*omega_n1;

w1 = (pi/L_tot)^2*sqrt(E*I/rho);
w2 = (2*pi/L_tot)^2*sqrt(E*I/rho);

% We now use V = PHI*S

% Load and free damping matrix
R_CDM = R0_static + R0_dynamic*(F1*sin(w1*t_CDM)+F2*sin(w2*t_CDM)); % Load vector som er N+1 bred
%R_CDM(:,3/dt+2:end) = 0;
Rf_CDM = R_CDM(df,:);
Rphi_CDM = PHi'*Rf_CDM;

Mff_phi = PHi'*Mff*PHi;
Cff_phi = zeros(size(Mff_phi,1),size(Mff_phi,2));
zeta = 0.20;
if c_damp > 0
  for i = 1:N_phi
    Cff_phi(i,i) = 2*zeta*m_phi(i,i)*Omega(i);
  end
end
Cff_phys = (Mff*PHi/Mff_phi)*Cff_phi*(Mff_phi\PHi'*Mff);

% Find the static displacement
[V0, ~, ~, ~, ~, ~] = NewmarkStatic_R0(X,T,D,G,U,R0_static,nno,nel,nd);

% setup initial conditions
%V0 = zeros(size(phi,1),1);
%S0 = zeros(size(PHi,2),1);
S0 = PHI\V0;
DV0 = zeros(size(phi,1),1);
DS0 = zeros(size(PHi,2),1);
Q_CDM = zeros(size(V0,1),N+1);
Q0 = Q(X,T,G,nel,nno,V0);
DDV0 = Mff\(Rf_CDM(:,1) - Cff_phys*DV0(df) - Q0(df));
DDS0 = (PHi'*Mff*PHi)\(Rphi_CDM(:,1) - Cff_phi*DS0 - PHi'*Q0(df));
V_before = V0(df) - dt*DV0(df) + (1/2)*dt^2*DDV0;
S_before = S0 - dt*DS0 + (1/2)*dt^2*DDS0; 


% Displacement matrix used in the CDM
V_CDM_Full = zeros(nd,N+1); % Matrix with the full degrees of freedom
V_CDM = [V_before, V0(df), zeros(size(V0(df),1),N-1)]; % matrix with only free degrees of freedom
V_CDM_Full(df,:) = V_CDM; % is the relation between the two matrices
Vs_CDM_Full = zeros(nd,N+1);
Vs_CDM = [V_before, V0(df), zeros(size(V0(df),1),N-1)];
Vs_CDM_Full(df,:) = Vs_CDM;
S_CDM = [S_before, S0, zeros(size(PHi,2),N-1)];

a_CDM_phi = (1/dt^2)*PHi'*Mff*PHi - (1/(2*dt))*Cff_phi;
b_CDM_phi = (2/dt^2)*PHi'*Mff*PHi;
Z_CDM_phi = (1/dt^2)*PHi'*Mff*PHi + (1/(2*dt))*Cff_phi;
a_CDM = (1/dt^2)*Mff - (1/(2*dt))*Cff_phys;
b_CDM = (2/dt^2)*Mff;
Z_CDM = (1/dt^2)*Mff + (1/(2*dt))*Cff_phys;

for i = 2:N 

  Qsi = Q(X,T,G,nel,nno,Vs_CDM_Full(:,i));
  Qsi = Qsi(df);
  Qi = Q(X,T,G,nel,nno,V_CDM_Full(:,i));
  Qi = Qi(df);

  S_CDM(:,i+1) = Z_CDM_phi\(Rphi_CDM(:,i-1) - PHi'*Qsi + b_CDM_phi*S_CDM(:,i) - a_CDM_phi*S_CDM(:,i-1));
  Vs_CDM_Full(df,i+1) = PHi*S_CDM(:,i+1);

  V_CDM(:,i+1) = Z_CDM\(Rf_CDM(:,i-1) - Qi + b_CDM*V_CDM(:,i) - a_CDM*V_CDM(:,i-1));
  V_CDM_Full(df,i+1) = V_CDM(:,i+1);

end


% testing the relationships between the modal derivatives
if N_phi == 3
  s4_test = S_CDM(2,:).*S_CDM(1,:);
  s5_test = (S_CDM(2,:)).^2 -(1/4)*(S_CDM(1,:)).^2;
  s6_test = S_CDM(3,:).*S_CDM(1,:) + (1/6)*(S_CDM(1,:)).^2;
  s7_test = S_CDM(3,:).*S_CDM(2,:);
  s8_test = (S_CDM(3,:)).^2;
end

%% Plot of midpoint using modal basis projection

unon_plot = V_CDM_Full(11*3-1,2:end);
us_plot = Vs_CDM_Full(11*3-1,2:end);

if V_plot == 1

  figure()
  plot(t_CDM,us_plot*1e2,'--k')
  hold on
  plot(t_CDM,unon_plot*1e2,':k.')
  xlabel('time, t [s]')
  ylabel('Displacement, u [cm]')
  title('Displacement over time using modal basis projection')
  legend({'CDM modal projection, u_{MD}','CDM nonlinear u_{non}'},'Location','northeast')
  hold off
end

if s_relations == 1
  figure()
  plot(t_CDM,S_CDM(4,2:end),':k.')
  hold on
  plot(t_CDM,s4_test(2:end),'--k')
  xlabel('time, t [s]')
  ylabel('modal coordinate s4')
  title('s4 vs s1^2-1/4*s2^2')
  legend({'s_4', 's_1^2-1/4*s_2^2 '},'Location','northeast')
  hold off

  figure()
  plot(t_CDM,S_CDM(5,2:end),':k.')
  hold on
  plot(t_CDM,s5_test(2:end),'--k')
  xlabel('time, t [s]')
  ylabel('modal coordinate s5')
  title('s5 vs s1*s2')
  legend({'s_5', 's_1*s_2 '},'Location','northeast')
  hold off

  figure()
  plot(t_CDM,S_CDM(6,2:end),':k.')
  hold on
  plot(t_CDM,s6_test(2:end),'--k')
  xlabel('time, t [s]')
  ylabel('modal coordinate s6')
  title('s6 vs s1*s3+3/8*s2^2')
  legend({'s_6', 's_1*s_3+3/8*s_2^2 '},'Location','northeast')
  hold off

  figure()
  plot(t_CDM,S_CDM(7,2:end),':k.')
  hold on
  plot(t_CDM,s7_test(2:end),'--k')
  xlabel('time, t [s]')
  ylabel('modal coordinate s7')
  title('s7 vs s2*s3')
  legend({'s_7', 's_2*s_3 '},'Location','northeast')
  hold off

  figure()
  plot(t_CDM,S_CDM(8,2:end),':k.')
  hold on
  plot(t_CDM,s8_test(2:end),'--k')
  xlabel('time, t [s]')
  ylabel('modal coordinate s8')
  title('s8 vs s3^2')
  legend({'s_8', 's_3^2 '},'Location','northeast')
  hold off
end

if Q_plot == 1
  figure()
  plot(V_CDM_Full(11*3-1,2:end),Q_CDM(11*3-1,2:end),':k.')
  hold on
  plot(V_CDM_Full_lin(11*3-1,2:end),Q_CDM_lin(11*3-1,2:end),'--k')
  xlabel('Displacement, u [m]')
  ylabel('Internal force, Q [N]')
  title('Internal force vs displacement')
  legend({'CDM non-linear', 'Q_{non},CDM linear, Q_{lin}'},'Location','northeast')
  hold off
end


%% functions

function Q1 = Q(X,T,G,nel,nno,V)
% Setup system force vector Q
  Q1 = zeros(nno*3,1);                   % initialization of Q1
  for el = 1:nel  
    no1 = T(el,1);  no2 = T(el,2);       % start node/end node
    X1 = X(no1,:);  X2 = X(no2,:);       % coordinates to start/end node
    de = [no1*3-2 no1*3-1 no1*3 ...
          no2*3-2 no2*3-1 no2*3];        % dofs for element
    ve = V(de);                          % displacements from system vector
    q = qbeam(X1,X2,G(el,:),ve);         % force vector for element
    Q1(de) = Q1(de) + q;                 % q is added to Q1
  end
end

function q = qbeam(X1,X2,Ge,ve)
% Setup element force vector
% integration points and weights

nip = 3; % Can be used up to polynomial order = 2*nip-1 = 5
zp = sqrt(3/5);
zip = [1-zp 1 1+zp]/2;
wip = [5/9 8/9 5/9]/2;

EA = Ge(1)*Ge(2);
EI = Ge(1)*Ge(3);
D = [EA 0; 0 EI];                      % material stiffness matrix
I1 = [1;0];
[A, L] = Abeam(X1,X2);                 % transformation matrix
vel = A*ve;                            % local displacement vector
q = zeros(6,1);

for ip=1:nip
  z = zip(ip);
  w = wip(ip);
  B = Bbeam(L,z);
  G = Gbeam(L,z);
  B1 = B+I1*(vel'*G);
  % strains
  e1 = B*vel+1/2*I1*(vel'*G*vel);
  % stresses 
  s = D*e1;
  % make q for local directions
  q = q+w*L*B1'*s;
end

% transform q to global directions
q = A'*q;
end

function K = K(X,T,G,nel,nno,V)
% Setup system stiffnessmatrix
K = zeros(nno*3,nno*3);                % initialization of K
for el = 1:nel
  no1 = T(el,1);  no2 = T(el,2);       % start node/end node
  X1 = X(no1,:);  X2 = X(no2,:);       % coordinates to start/end node
  de = [no1*3-2 no1*3-1 no1*3 ...
        no2*3-2 no2*3-1 no2*3];        % dofs for element
  ve = V(de);                          % displacements from system vector
  k = kbeam(X1,X2,G(el,:),ve);         % stiffnessmatrix for element
  K(de,de) = K(de,de) + k;             % k is added to K
end
end

function k = kbeam(X1,X2,Ge,ve)
% Setup non-linear element stiffness matrix 

% integration points and weights
nip = 3;
zp = sqrt(3/5);
zip = [1-zp 1 1+zp]/2;
wip = [5/9 8/9 5/9]/2;

EA = Ge(1)*Ge(2);
EI = Ge(1)*Ge(3);
D = [EA 0; 0 EI];                      % material stiffness matrix
I1 = [1;0];
[A, L] = Abeam(X1,X2);                  % transformation matrix
vel = A*ve;                            % local displacement vector
k = zeros(6);
for ip=1:nip
  z = zip(ip);
  w = wip(ip);
  B = Bbeam(L,z);
  G = Gbeam(L,z);
  % strain
  e1 = B*vel+1/2*I1*(vel'*G*vel);
  % stress (normalforce)
  s = EA*e1(1);
  % make k for local directions
  k = k+w*L*s*G;
  B1 = B+(1/2)*I1*vel'*G;
  k = k+w*L*B'*D*B1;
end
% transform k to global directions
k = A'*k*A;
end

function Kt = Kt(X,T,G,nel,nno,V)
% Setup system stiffnessmatrix
Kt = zeros(nno*3,nno*3);                % initialization of K
for el = 1:nel
  no1 = T(el,1);  no2 = T(el,2);       % start node/end node
  X1 = X(no1,:);  X2 = X(no2,:);       % coordinates to start/end node
  de = [no1*3-2 no1*3-1 no1*3 ...
        no2*3-2 no2*3-1 no2*3];        % dofs for element
  ve = V(de);                          % displacements from system vector
  kt = ktbeam(X1,X2,G(el,:),ve);         % stiffnessmatrix for element
  Kt(de,de) = Kt(de,de) + kt;             % k is added to K
end
end

function kt = ktbeam(X1,X2,Ge,ve)
% Setup element stiffness matrix - tangent

% integration points and weights
nip = 3;
zp = sqrt(3/5);
zip = [1-zp 1 1+zp]/2;
wip = [5/9 8/9 5/9]/2;

EA = Ge(1)*Ge(2);
EI = Ge(1)*Ge(3);
D = [EA 0; 0 EI];                      % material stiffness matrix
I1 = [1;0];
[A, L] = Abeam(X1,X2);                  % transformation matrix
vel = A*ve;                            % local displacement vector
kt = zeros(6);
for ip=1:nip
  z = zip(ip);
  w = wip(ip);
  B = Bbeam(L,z);
  G = Gbeam(L,z);
  % strain
  e1 = B*vel+1/2*I1*(vel'*G*vel);
  % stress (normalforce)
  s = EA*e1(1);
  % make k for local directions
  kt = kt+w*L*s*G;
  B1 = B+I1*(vel'*G);
  kt = kt+w*L*B1'*D*B1;
end
% transform k to global directions
kt = A'*kt*A;
end

function dKs = dKs(X,T,G,nel,nno,phi,V)
% Setup system stiffnessmatrix
dKs = zeros(nno*3,nno*3);                % initialization of K
for el = 1:nel
  no1 = T(el,1);  no2 = T(el,2);       % start node/end node
  X1 = X(no1,:);  X2 = X(no2,:);       % coordinates to start/end node
  de = [no1*3-2 no1*3-1 no1*3 ...
        no2*3-2 no2*3-1 no2*3];        % dofs for element
  ve = V(de);                          % displacements from system vector
  phie = phi(de);
  dks = dksbeam(X1,X2,G(el,:),phie,ve);         % stiffnessmatrix for element
  dKs(de,de) = dKs(de,de) + dks;             % k is added to K
end
end

function dks = dksbeam(X1,X2,Ge,phie,ve)
% Setup element stiffness matrix - tangent

% integration points and weights
nip = 3;
zp = sqrt(3/5);
zip = [1-zp 1 1+zp]/2;
wip = [5/9 8/9 5/9]/2;

EA = Ge(1)*Ge(2);
EI = Ge(1)*Ge(3);
[A, L] = Abeam(X1,X2);                  % transformation matrix
vel = A*ve;                     % local displacement vector
phiel = A*phie;

dks = zeros(6);
for ip=1:nip
  z = zip(ip);
  w = wip(ip);
  B = Bbeam(L,z);
  G = Gbeam(L,z);
  % B-matrix strain membrane term
  B_eps = B(1,:);
  % First term
  ef1 = B_eps'*phiel'*G;
  % Second term
  ef2 = 2*G*phiel*B_eps;
  % Third term
  ef3 = G*(phiel*vel' + vel*phiel')*G;
  % make k for local directions
  dks = dks+w*L*EA/2*(ef1 + ef2 + ef3);
end
% transform k to global directions
dks = A'*dks*A;
end

% Set up load vector R
function R = R(X,T,D,dL,bL,nd,nel)
R = zeros(nd,1);                       % initialization of R
for el = 1:nel
  dLe = dL(el,:);
  if dLe(1) > 0
    no1 = T(el,1);  no2 = T(el,2);     % startnode/endnode
    X1 = X(no1,:);  X2 = X(no2,:);     % coordinats to start-/endnode
    r = rbeam(X1,X2,dLe);              % element load vector
    de = D(el,:);                      % indexarray for element nodes
    R(de) = R(de) + r;                 % r placed into R
  end
end
nbL = size(bL,1);                      % number  of nodal loads
for i = 1:nbL
  d = bL(i,1);                         % global dof with nodal laod
  R(d) = R(d) + bL(i,2);               % nodal load placed into R 
end
end

function r = rbeam(X1,X2,dLe)
% Set up element load vector
% transformation matrix
[A, L] = Abeam(X1,X2);
% r in local system
r = zeros(6,1);
if dLe(1) == 1
  p = dLe(2)*L/2;
  r = [p 0 0 p 0 0]';
elseif dLe(1) == 2
  p = dLe(2)*L/2;
  m = dLe(2)*L^2/12;
  r = [0 p m 0 p -m]';
else
  disp('Fejl i specifikation af last!')
end
% transform r to global system
r = A'*r;
end

function M = M(X,T,G,D,nd,nel)
M = zeros(nd,nd);
for el = 1:nel
  no1 = T(el,1);  no2 = T(el,2);       % startnode/endnode
  X1 = X(no1,:);  X2 = X(no2,:);       % coordinats to start-/endnode
  m = mbeam(X1,X2,G(el,:));            % element stiffness matrix
  de = D(el,:);                        % indexarray for element nodes
  M(de,de) = M(de,de) + m;             % k placed into K
end
end

function m = mbeam(X1, X2,Ge)
%Transformation matrix og længde
[A, L] = Abeam(X1,X2);
%M i lokal system
rho = Ge(4);
m = rho*[L/3     0          0       L/6      0             0
         0    13*L/35   11*L^2/210   0     9*L/70      -13*L^2/420
         0   11*L^2/210  L^3/105     0   13*L^2/420    -L^3/140
         L/6     0          0       L/3      0             0
         0    9*L/70   13*L^2/420    0    13*L/35      -11*L^2/210
         0 -13*L^2/420  -L^3/140     0   -11*L^2/210       L^3/105];
%Transformer til gloabalt system
m = A'*m*A;
end

function C = C(X,T,G,D,nd,nel)
C = zeros(nd,nd);
for el = 1:nel
  no1 = T(el,1);  no2 = T(el,2);       % startnode/endnode
  X1 = X(no1,:);  X2 = X(no2,:);       % coordinats to start-/endnode
  c = cbeam(X1,X2,G(el,:));            % element stiffness matrix
  de = D(el,:);                        % indexarray for element nodes
  C(de,de) = C(de,de) + c;             % k placed into K
end
end

function c = cbeam(X1, X2, Ge)
%Transformation matrix og længde
[A, L] = Abeam(X1,X2);
%M i lokal system
c_damp = Ge(5);
c = c_damp*[L/3     0          0       L/6      0             0
         0    13*L/35   11*L^2/210   0     9*L/70      -13*L^2/420
         0   11*L^2/210  L^3/105     0   13*L^2/420    -L^3/140
         L/6     0          0       L/3      0             0
         0    9*L/70   13*L^2/420    0    13*L/35      -11*L^2/210
         0 -13*L^2/420  -L^3/140     0   -11*L^2/210       L^3/105];
%Transformer til gloabalt system
c = A'*c*A;
end

function G1 = Gbeam(L,z)
% setup strain interpolation matrix related to nonlinearity
  Gx = [-1/L 0 0 1/L 0 0];
  Gy = zeros(1,6);
  Gy(1,2) = -6*(z-z^2)/L;
  Gy(1,3) = 1-4*z+3*z^2;
  Gy(1,5) = 6*(z-z^2)/L;
  Gy(1,6) = -2*z+3*z^2;
  G1 = 0*(Gx'*Gx) + (Gy'*Gy); % second order x-term taken out
end

function B = Bbeam(L,z)
% setup standard linear strain interpolation matrix
  B = zeros(2,6);
  B(1,1) = -1/L;
  B(1,4) = 1/L;
  B(2,2) = -(6-12*z)/L^2;
  B(2,3) = (-4+6*z)/L;
  B(2,5) = (6-12*z)/L^2;
  B(2,6) = (-2+6*z)/L; 
end


function [A, L] = Abeam(X1,X2)
% Setup transformation matrix A
n = X2-X1;                      % direction vector
L = sqrt(dot(n,n));             % element length
n = n/L;                        % unit vector
% establish transformationsmatrix
A = [ n(1) n(2)   0     0    0    0
     -n(2) n(1)   0     0    0    0
       0     0    1     0    0    0
       0     0    0   n(1) n(2)   0
       0     0    0  -n(2) n(1)   0
       0     0    0     0    0    1];
end

function plotTop(X,T,nel,nno)
% plot geometry
figure(); title('Elementtopologi'); axis equal; hold on
for el = 1:nel
  plot(X(T(el,:),1),X(T(el,:),2),'b-')
end
for no = 1:nno
  text(X(no,1),X(no,2),num2str(no),...
       'color','blue','BackgroundColor',[0.7 0.7 0.7]);
end
for el = 1:nel
  xp=mean(X(T(el,:),:));
  text(xp(1),xp(2),num2str(el),'color','black');
end
hold off
end

function plotDof(X,T,D,V,nel,skala)
% plot deformed geometry (degrees of freedom)
figure(); axis equal; hold on 
title(['Deformationer, skala: ' num2str(skala, '%10.3e')]); 
for el = 1:nel
  plot(X(T(el,:),1),X(T(el,:),2),'b:')
  % deformed
  % set up transformation matrix
  no1 = T(el,1);  no2 = T(el,2);       % startknude/slutknude
  X1 = X(no1,:);  X2 = X(no2,:);       % koordinater til start-/slutknude
  [A L] = Abeam(X1,X2);
  % set up transformation matrix for displacements
  Au=A(1:2,1:2);
  % get lokal dofs
  v=V(D(el,:));
  % coordinates plus displacements
  nrp=11;
  Xs=zeros(2,nrp);
  for i=1:nrp
    s=(i-1)/(nrp-1);
    N=[1-s             0               0 s           0            0;
       0   1-3*s^2+2*s^3 (s-2*s^2+s^3)*L 0 3*s^2-2*s^3 (-s^2+s^3)*L];
    Xs(:,i)=X(T(el,1),:)'*(1-s)+X(T(el,2),:)'*s+skala*Au'*N*A*v;
  end
  plot(Xs(1,:),Xs(2,:),'b-');
end
hold off
end