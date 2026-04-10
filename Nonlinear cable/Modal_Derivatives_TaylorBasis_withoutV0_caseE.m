% Non-linear FEM (Modal derivatives Taylor - Lagrange) 

clear all
close all

%% Input data

% Geometric element parameters
d = 50e-3; % [m] diameter of the cable
A = (pi/4)*(d)^2; % [m^2] Cross sectional area
I = pi*(d)^4/64; % [m^4] second moment of inertia 
a_L = 1; % [m] length of element
L_tot = 20; % [m] Total length of element (use if easier)

% Material parameters
E = 210e9; % [Pa] Young's modulus
rho = 15.6; % [kg/m] Density of the element
c_damp = 0; % damping coefficient of the element

% Load parameters
F1 = -1; % [N] Point load
F2 = -3; % [N/m] Line load
e = 0; % [m] Eccentricity
w1 = (pi/L_tot)^2*sqrt(E*I/rho);
w2 = (2*pi/L_tot)^2*sqrt(E*I/rho);
p = rho*9.81;
P = 4.4237e+03;

% Prestressing
N0 = 4.4237e+03;

% Mode parameters
N_phi = 20; % Number of linear modes to use

% Damping ratio
zeta = 0.20;

% Plotting parameters
Geo_plot = 0;

% Switch for using MD's
MD_phi = 1;

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

G = zeros(nel,6);

for i = 1:nel
  G(i,:) = [E A I rho c_damp N0];
end

% understøtninger: U(i)=global dof
U(1) = 1;
U(2) = 2;
U(3) = nno*3-1;
%U(4) = nno*3-1;

% nodal load (boundary Load): bL(i,:)=[global_dof, magnitude]
%bL = []; %empty if no load
bL = [nno*3-2 P];
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
if Geo_plot == 1
  plotTop(X,T,nel,nno)
end

%% Program

% setup constant coefficient matrices of the system
M_mass = M(X,T,G,D,nd,nel); % Mass matrix
C_damp = C(X,T,G,D,nd,nel); % Damping matrix
R0_dynamic = R(X,T,D,dL_dynamic,[],nd,nel); % Load amplitude vector
R0_static = R(X,T,D,dL_static,bL,nd,nel);
Q0 = zeros(nd,1);
q0 = [-1 0 0 1 0 0]';
for el = 1:nel
  N0 = G(el,6);
  q0l = N0*q0;
  no1 = T(el,1);  no2 = T(el,2);     
  X1 = X(no1,:);  X2 = X(no2,:); 
  [A, L] = Abeam(X1,X2);
  de = D(el,:);
  Q0(de) = Q0(de) + A'*q0l;
end
R_static = R0_static - Q0;

% setup displacement vector
V = zeros(nno*3,1);

K0 = Kt(X,T,G,nel,nno,V);


%% solve linear eigenvalue problem for the EOM
dof = 1:3*nno;                         % index to all dofs 
du = U;                                % index to prescribed dofs
df = setdiff(dof,du);                  % index to free dofs

K0ff = K0(df,df);                          % Free dofs only 
Mff = M_mass(df,df);                          % Free dofs only

% Modes and frequencies for zero initial displacement

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


%% Solving the complete modal derivatives equation

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


% Solving the complete modal derivative equation for proper damping of the
% full physical response simulation

HH = size(phi,2)*(size(phi,2)+1)/2;
dPhi_phys = zeros(nd,HH);

for i = 1:size(phi,2)
  for j = 1:i
    dKsi = dKs(X,T,G,nel,nno,phi(:,i),V);
    dKsj = dKs(X,T,G,nel,nno,phi(:,j),V);
    alpha0 = phi(df,j)'/m_phi(j,j)*(dKsi(df,df)*phi(df,j) - dKsj(df,df)*phi(df,i));
    dphi_phys = (K0ff - (Omega(i)+Omega(j))^2*Mff)\((alpha0*Mff - dKsi(df,df))*phi(df,j)*2);
    if i==j
      dphi_phys = dphi_phys/2;
    end

    % Map (i,j) to global modal derivative index
    k = i*(i-1)/2 + j;

    dPhi_phys(df,k) = dphi_phys;  
  end
end

%% Setting up the physical modal basis for damping

% Initialize the PHI-matrix containing the modes and the modal derivatives

PHI_phys = [phi dPhi_phys];
PHi_phys = [phi(df,:) dPhi_phys(df,:)];

orig_idx = 1:size(PHI_phys,2);
dependent_idx = [];

n0 = size(PHI_phys,2);
Beta = NaN(n0, n0);    % now use n0 x n0 for original-coordinate storage

if MD_phi == 1
    i = size(phi,2) + 1;

    while i <= size(PHI_phys, 2)
        PHI_test = PHI_phys;
        PHI_test(:,i) = [];

        %beta = (PHI_test' * PHI_test) \ (PHI_test' * PHI_phys(:,i));
        beta = lsqminnorm(PHI_test, PHI_phys(:,i));
        Top = PHI_phys(:,i) - PHI_test * beta;
        e_error = norm(Top) / norm(PHI_phys(:,i)) * 100;

        if e_error <= 1
            dep_col = orig_idx(i);
            dependent_idx(end+1) = dep_col;

            % original column indices still present in PHI_test
            remaining_idx = orig_idx;
            remaining_idx(i) = [];

            beta_full = NaN(n0,1);
            beta_full(remaining_idx) = beta;

            Beta(:,dep_col) = beta_full;

            PHI_phys(:,i) = [];
            orig_idx(i) = [];
        else
            i = i + 1;
        end
    end

    PHi_phys = PHI_phys(df,:);
else
    PHI_phys = phi;
    PHi_phys = phi(df,:);
end


%% setting up the Taylor basis matrix
% Initialize the PHI-matrix containing the modes and the modal derivatives


PHI = [phi(:,1:N_phi) dPhi(:,1:end)]; % Phi-matrix with zeros for the nonfree dofs
PHi = [phi(df,1:N_phi) dPhi(df,1:end)];% phi-matrix only containing the free dofs

%% Set-up of iteration parameters

dt = 1e-4;
t_CDM = 0:dt:30;
N = length(t_CDM);

%% solve for independent s-functions using CDM

% We now use V = PHI*S

% Load vector
t_shift = 20;
R_dynamic = R0_dynamic*(F1*sin(w1*(t_CDM-t_shift))+F2*sin(w2*(t_CDM-t_shift)));
R_dynamic(:,t_CDM<t_shift) = 0;
R_CDM = (R0_static+R_dynamic) - Q0; % Load vector som er N+1 bred
Rf_CDM = R_CDM(df,:);

% Find the static displacement
%[V0, ~, ~, ~, ~, ~, ~, ~] = NewmarkStatic_taylor(X,T,D,G,U,bL,dL,nno,nel,nd,PHI,N_phi);
%ST0 = Qbeam(N_phi,s0)*s0;

% setup initial conditions
V0 = zeros(size(PHI,1),1); %Physical initial displacement of the system
DV0 = zeros(size(PHI,1),1); %Physical initial velocity of the system

% Initial taylor modal displacement
ST0 = zeros(size(PHi,2),1);


% New coordinate s initial states
s0 = zeros(N_phi,1);
ds0 = zeros(N_phi,1);

% Set up modal matrices and external loading for taylor basis
mT = PHi'*Mff*PHi;
cT_red = zeros(size(mT,1),size(mT,2)); % damping on linear modes = damping on MDs (setup)

if zeta > 0
  for i = 1:N_phi
    cT_red(i,i) = 2*zeta*m_phi(i,i)*Omega(i);
  end
end
fT = PHi'*Rf_CDM;

% Set up initial internal force
Qs0 = Q(X,T,G,nel,nno,zeros(nd,1));
Qs0 = Qs0(df);
ghat0 = PHi'*Qs0;
Us0 = Ubeam(N_phi,s0);
gbar0 = Us0'*ghat0;

% Evaluate initial acceleration for new coordinate s
Pds0 = Pbeam(N_phi,ds0);
ms0 = Us0'*mT*Us0;
csds0 = Us0'*(2*mT*Pds0 + cT_red*Us0);
fst0 = Us0'*fT(:,1);
dds0 = ms0\(fst0 - csds0*ds0 - gbar0);

% Evaluate modal displacements from taylor series
s_before1 = s0 + (-dt)*ds0 + (1/2)*(-dt)^2*dds0;
s_before2 = s0 + (-2*dt)*ds0 + (1/2)*(-2*dt)^2*dds0;

% Set up storage of calculated values
VsT_CDM_Full = zeros(nd,N);
VsT_CDM = [V0(df), zeros(size(V0(df),1),N-1)];
VsT_CDM_Full(df,:) = VsT_CDM;
ST_CDM = [ST0, zeros(size(PHi,2),N-1)];
s_CDM = [s_before1, s0, zeros(N_phi,N-1)];
ds_CDM = zeros(size(s_CDM,1),size(s_CDM,2));
dST_CDM = zeros(size(ST_CDM,1),size(ST_CDM,2));
dVsT_CDM_Full = zeros(nd,N);

for i = 2:N
  %Evaluate approximated modal velocities
  if i == 2
    ds_CDM(:,i) = (1/dt)*(3/2*s_CDM(:,i) - 2*s_CDM(:,i-1) + 1/2*s_before2);
  else
    ds_CDM(:,i) = (1/dt)*(3/2*s_CDM(:,i) - 2*s_CDM(:,i-1) + 1/2*s_CDM(:,i-2));
  end
  %Evaluate transformation matrices
  Us = Ubeam(N_phi,s_CDM(:,i));
  Pds = Pbeam(N_phi,ds_CDM(:,i));
  %Evaluate velocity
  dST_CDM(:,i) = Us*ds_CDM(:,i);
  %Evaluate mass, damping and external loading
  ms = Us'*mT*Us;
  csds = Us'*(2*mT*Pds + cT_red*Us);
  fst = Us'*fT(:,i-1);
  %Evaluate internal restoring forces
  Qsi = Q(X,T,G,nel,nno,VsT_CDM_Full(:,i-1));
  Qsi = Qsi(df);
  ghat = PHi'*Qsi;
  gbar = Us'*ghat;
  %Evaluate next modal displacement by CDM
  ZT = ms/dt^2 + csds/(2*dt);
  bT = 2/dt^2*ms;
  aT = ms/dt^2 - csds/(2*dt);
  s_CDM(:,i+1) = ZT\(fst - gbar + bT*s_CDM(:,i) - aT*s_CDM(:,i-1));
  %Transform modal solution into physical coordinates
  Qs = Qbeam(N_phi,s_CDM(:,i+1));
  ST_CDM(:,i) = Qs*s_CDM(:,i+1);
  VsT_CDM_Full(df,i) = PHi*ST_CDM(:,i) + V0(df);
  dVsT_CDM_Full(df,i) = PHi*dST_CDM(:,i); 
end

V_phys = V_physical_pre(X,T,G,D,R_CDM,U,nno,nel,nd,dt,t_CDM,V0,DV0,zeta,PHI_phys,size(phi,2),Omega);


%% Plots
u_nonst = VsT_CDM_Full(11*3-1,:); %Taylor nonlinear displacement at midpoint by taylor basis
u_nonlinear = V_phys(11*3-1,2:end);

% Taylor basis
figure()
plot(t_CDM,u_nonst*1e2,'-k')
hold on
plot(t_CDM,u_nonlinear*1e2,':k.')
xlabel('Time, t [s]')
ylabel('displacement, u [cm]')
title('taylor basis displacement of midpoint')
legend({'CDM nonlinear taylor basis, u_{nonsT}','CDM nonlinear u_{non}'},'Location','northeast')


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
N0 = Ge(6);
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
  s = D*e1; %+ [N0 0]';
  % make q for local directions
  q = q+w*L*(B1'*s + N0*G*vel);
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
N0 = Ge(6);
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
  s = EA*e1(1) + N0;
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
N0 = Ge(6);
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
  s = EA*e1(1)+N0;
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

function dKs_V0 = dKs_V0(X,T,G,nel,nno,phi,V)
% Setup system stiffnessmatrix
dKs_V0 = zeros(nno*3,nno*3);                % initialization of K
for el = 1:nel
  no1 = T(el,1);  no2 = T(el,2);       % start node/end node
  X1 = X(no1,:);  X2 = X(no2,:);       % coordinates to start/end node
  de = [no1*3-2 no1*3-1 no1*3 ...
        no2*3-2 no2*3-1 no2*3];        % dofs for element
  ve = V(de);                          % displacements from system vector
  phie = phi(de);
  dks_V0 = dksbeam_V0(X1,X2,G(el,:),phie,ve);         % stiffnessmatrix for element
  dKs_V0(de,de) = dKs_V0(de,de) + dks_V0;             % k is added to K
end
end

function dks_V0 = dksbeam_V0(X1,X2,Ge,phie,ve)
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

dks_V0 = zeros(6);
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
  ef3 = G*(2*phiel*vel' + vel*phiel')*G;
  % make k for local directions
  dks_V0 = dks_V0+w*L*EA/2*(ef1 + ef2 + ef3);
end
% transform k to global directions
dks_V0 = A'*dks_V0*A;
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
  G1 = (Gx'*Gx) + (Gy'*Gy); % second order x-term taken out
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

function Us = Ubeam(N,s)

H = Hbeam(N);
Bs = zeros(size(H,1),size(H,2));

for k = 1:N
  Bk = Bkbeam(N,k);
  Bs = Bs + Bk*s(k);
end

Us = H + Bs;

end

function Pds = Pbeam(N,ds)
H = N*(N+1)/2;
Pds = zeros(N+H,N);
for k = 1:N
  Ak = Akbeam(N,k);
  Pds = Pds + Ak*ds(k);
end

end

function Qs = Qbeam(N,s)
H = Hbeam(N);
Ws = zeros(size(H,1),size(H,2));
for k = 1:N
  Ak = Akbeam(N,k);
  Ws = Ws + Ak*s(k);
end
Qs = H + Ws;
end

function H = Hbeam(N)
% Computes the matrix H for taylor basis MD's
% Compute the relevant identity matrix
I = eye(N);

% compute the number of MD's for the system
HH = N*(N+1)/2;

% make the H matrix
H = zeros(N+HH,N);
int = 1:N;
H(int,int) = I;

end

function Bk = Bkbeam(N, k)

    % --- Dimensions ---
    H  = N*(N+1)/2;
    Bk = Akbeam(N, k);

    % --- Build mapping from (p,q) to row index in the d-vector ---
    % Ordering of pairs (p,q): (1,1), (2,1),(2,2), (3,1),(3,2),(3,3), ...
    pairRow = zeros(N,N);
    row = 1;
    for p = 1:N
        for q = 1:p
            pairRow(p,q) = row;
            row = row + 1;
        end
    end

    % --- Add contributions of e_{ik} v_i^T ---
    %
    % j = k is fixed.
    %
    for i = k:N
        r = N + pairRow(i, k);   % row index inside the (N+H)-vector
        Bk(r, i) = Bk(r, i) + 1; % add the e_{ik} v_i^T term
    end
end

function A_k = Akbeam(N,k)
% A_k is (N+H)-by-N   where H = N(N+1)/2.

H = N*(N+1)/2;
A_k = zeros(N+H, N);

% Triangular index base offset
base = (k-1)*k/2;

% Positions of d_{k,1}, d_{k,2}, ..., d_{k,k} in the last H entries
idx = N + (base + (1:k));

% Fill A_k row blocks:
% A_k(idx(j), j) = 1 for j = 1..k
A_k(sub2ind([N+H, N], idx, 1:k)) = 1;

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