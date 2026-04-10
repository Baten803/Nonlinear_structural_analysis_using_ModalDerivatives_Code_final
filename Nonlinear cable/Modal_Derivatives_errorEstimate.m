% Non-linear FEM (Modal derivatives - Lagrange) 
tic
clear all
close all

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

% Switches
linear_phi = 1;
MD = 0;
V_plot = 1;
Q_plot = 0;
s_relations = 1;


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
 dL = zeros(nel,2);
for i = 11:nel
dL(i,:) = [2 -1];
end 

%% Program

% setup displacement vector
V = zeros(nno*3,1);

% setup constant coefficient matrices of the system
M_mass = M(X,T,G,D,nd,nel); % Mass matrix
C_damp = C(X,T,G,D,nd,nel); % Damping matrix
K0 = K(X,T,G,nel,nno,V); % standard linear stiffness matrix
R0 = R(X,T,D,dL,bL,nd,nel); % Load amplitude vector

%% Set-up of iteration parameters
dof = 1:3*nno;                         % index to all dofs 
du = U;                                % index to prescribed dofs
df = setdiff(dof,du);                  % index to free dofs

T_Program = toc;

tic
dt_phys = 1e-4;
t_CDM = 0:dt_phys:10;
N = length(t_CDM);

%% solve for dynamic response V

w1 = (pi/L_tot)^2*sqrt(E*I/rho);
w2 = (2*pi/L_tot)^2*sqrt(E*I/rho);

% We now use V = PHI*S

% Load and free degree matrices 
R_CDM = R0*(F1*sin(w1*t_CDM)+F2*sin(w2*t_CDM)); % Load vector som er N+1 bred
Rf_CDM = R_CDM(df,:);
Mff = M_mass(df,df);
Cff = C_damp(df,df);

% setup initial conditions
V0 = zeros(nd,1);
DV0 = zeros(nd,1);
Q0 = Q(X,T,G,nel,nno,V0);
DDV0 = Mff\(Rf_CDM(:,1) - Cff*DV0(df) - Q0(df));
V_before = V0(df) - dt_phys*DV0(df) + (1/2)*dt_phys^2*DDV0; 

% Displacement matrix used in the CDM
V_CDM_Full = zeros(nd,N+1); % Matrix with the full degrees of freedom
V_CDM = [V_before, V0(df), zeros(size(V0(df),1),N-1)]; % matrix with only free degrees of freedom
V_CDM_Full(df,:) = V_CDM; % is the relation between the two matrices


a_CDM = (1/dt_phys^2)*Mff - (1/(2*dt_phys))*Cff;
b_CDM = (2/dt_phys^2)*Mff;
Z_CDM = (1/dt_phys^2)*Mff + (1/(2*dt_phys))*Cff;

for i = 2:N 

  Qi = Q(X,T,G,nel,nno,V_CDM_Full(:,i));
  Qi = Qi(df);

  V_CDM(:,i+1) = Z_CDM\(Rf_CDM(:,i-1) - Qi + b_CDM*V_CDM(:,i) - a_CDM*V_CDM(:,i-1));
  V_CDM_Full(df,i+1) = V_CDM(:,i+1);

end

T_V_physical_ref = toc;

%% Physical coordinates

dt_phys_lim = 1.13e-4;
t_CDM_physlim = 0:dt_phys_lim:10;
R_CDM_physlim = R0*(F1*sin(w1*t_CDM_physlim)+F2*sin(w2*t_CDM_physlim));

tic
V_phys_lim = V_physical(X,T,G,D,R_CDM_physlim,U,nno,nel,nd,dt_phys_lim,V0,DV0);
T_V_physical_lim = toc;

%% Standard basis approximation

MD_no = 0;

% 2 modes
Ns_phi_2 = 2;
Hs_phi_2 = Ns_phi_2*(Ns_phi_2+1)/2;
S0_2 = zeros(Ns_phi_2+Hs_phi_2,1);
DS0_2 = zeros(Ns_phi_2+Hs_phi_2,1);
dt_s2 = dt_phys;%6.09e-4;
t_CDM_s2 = 0:dt_s2:10;
R_CDM_s2 = R0*(F1*sin(w1*t_CDM_s2)+F2*sin(w2*t_CDM_s2));

% 3 modes
Ns_phi_3 = 3;
Hs_phi_3 = Ns_phi_3*(Ns_phi_3+1)/2;
S0_3 = zeros(Ns_phi_3+Hs_phi_3,1);
DS0_3 = zeros(Ns_phi_3+Hs_phi_3,1);
dt_s3 = dt_phys;%3.979e-4;
t_CDM_s3 = 0:dt_s3:10;
R_CDM_s3 = R0*(F1*sin(w1*t_CDM_s3)+F2*sin(w2*t_CDM_s3));

tic
Vsb_2_3 = V_StandardBasis(X,T,G,D,R_CDM_s2,U,nno,nel,nd,Ns_phi_2,MD_no,dt_s2,S0_2,DS0_2);
T_Vsb_2_3 = toc;

tic
Vsb_3_5 = V_StandardBasis(X,T,G,D,R_CDM_s3,U,nno,nel,nd,Ns_phi_3,MD_no,dt_s3,S0_3,DS0_3);
T_Vsb_3_5 = toc;

%% Taylor basis approximation

c_damp = 0;

% 2 modes
Nt_phi_2 = 2;
Ht_phi_2 = Nt_phi_2*(Nt_phi_2+1)/2;
ST0_2 = zeros(Nt_phi_2+Ht_phi_2,1);
DST0_2 = zeros(Nt_phi_2+Ht_phi_2,1);
dt_t2 = 10e-3;%0.179; 
t_CDM_t2 = 0:dt_t2:10;
R_CDM_t2 = R0*(F1*sin(w1*t_CDM_t2)+F2*sin(w2*t_CDM_t2));

% 3 modes
Nt_phi_3 = 3;
Ht_phi_3 = Nt_phi_3*(Nt_phi_3+1)/2;
ST0_3 = zeros(Nt_phi_3+Ht_phi_3,1);
DST0_3 = zeros(Nt_phi_3+Ht_phi_3,1);
dt_t3 = 10e-3;%0.13;
t_CDM_t3 = 0:dt_t3:10;
R_CDM_t3 = R0*(F1*sin(w1*t_CDM_t3)+F2*sin(w2*t_CDM_t3));

tic
Vtb_2 = V_TaylorBasis(X,T,G,D,R_CDM_t2,U,nno,nel,nd,Nt_phi_2,c_damp,dt_t2,ST0_2,DST0_2);
T_Vtb_2 = toc;

tic
Vtb_3 = V_TaylorBasis(X,T,G,D,R_CDM_t3,U,nno,nel,nd,Nt_phi_3,c_damp,dt_t3,ST0_3,DST0_3);
T_Vtb_3 = toc;

%% Time for running the functions
fprintf('Program :%.6f s\n', T_Program');
fprintf('V_physical    : %.6f s\n', T_V_physical_ref);
fprintf('V_physical    : %.6f s\n', T_V_physical_lim);
fprintf('V_StandardBasis (2 modes): %.6f s\n', T_Vsb_2_3);
fprintf('V_StandardBasis (3 modes): %.6f s\n', T_Vsb_3_5);
fprintf('V_TaylorBasis   (2 modes): %.6f s\n', T_Vtb_2);
fprintf('V_TaylorBasis   (3 modes): %.6f s\n', T_Vtb_3);

%% Error estimation

%Physical reference response at midpoint
V_phys_mid = V_CDM_Full(11*3-1,2:end);

%Physical limit response at midpoint
V_phys_mid_lim = V_phys_lim(11*3-1,2:end);

% Only interpolate on the overlapping time interval
idx_phys_lim = (t_CDM >= t_CDM_physlim(1)) & (t_CDM <= t_CDM_physlim(end));

V_phys_mid_lim_interp = nan(size(t_CDM));
V_phys_mid_lim_interp(idx_phys_lim) = interp1(t_CDM_physlim, V_phys_mid_lim, t_CDM(idx_phys_lim), 'linear');

%Standard basis responses at midpoint
Vsb_2_3_mid = Vsb_2_3(11*3-1,2:end);
Vsb_3_5_mid = Vsb_3_5(11*3-1,2:end);

% Only interpolate on overlapping time interval
idx_s2 = (t_CDM >= t_CDM_s2(1)) & (t_CDM <= t_CDM_s2(end));
idx_s3 = (t_CDM >= t_CDM_s3(1)) & (t_CDM <= t_CDM_s3(end));

Vsb_2_mid_interp = nan(size(t_CDM));
Vsb_3_mid_interp = nan(size(t_CDM));

Vsb_2_mid_interp(idx_s2) = interp1(t_CDM_s2, Vsb_2_3_mid, t_CDM(idx_s2), 'linear');
Vsb_3_mid_interp(idx_s3) = interp1(t_CDM_s3, Vsb_3_5_mid, t_CDM(idx_s3), 'linear');


%Taylor basis responses at midpoint
Vtb_2_mid = Vtb_2(11*3-1,:);
Vtb_3_mid = Vtb_3(11*3-1,:);

% Only interpolate on overlapping time interval
idx_t2 = (t_CDM >= t_CDM_t2(1)) & (t_CDM <= t_CDM_t2(end));
idx_t3 = (t_CDM >= t_CDM_t3(1)) & (t_CDM <= t_CDM_t3(end));

Vtb_2_mid_interp = nan(size(t_CDM));
Vtb_3_mid_interp = nan(size(t_CDM));

Vtb_2_mid_interp(idx_t2) = interp1(t_CDM_t2, Vtb_2_mid, t_CDM(idx_t2), 'linear');
Vtb_3_mid_interp(idx_t3) = interp1(t_CDM_t3, Vtb_3_mid, t_CDM(idx_t3), 'linear');

%Error estimate for physical coordinate limit
Vdiff_phys_lim = V_phys_mid(idx_phys_lim) - V_phys_mid_lim_interp(idx_phys_lim);

e_phys_lim_MAX = max(abs(Vdiff_phys_lim)/max(abs(V_phys_mid(idx_phys_lim))));
e_phys_lim_AVG = (sum(abs(Vdiff_phys_lim))/length(Vdiff_phys_lim))/(sum(abs(V_phys_mid(idx_phys_lim)))/length(V_phys_mid(idx_phys_lim)));

%Error estimate for standard basis
V_phys_mid_s2 = V_phys_mid(idx_s2);
V_phys_mid_s3 = V_phys_mid(idx_s3);

Vdiff_standard_2 = V_phys_mid_s2 - Vsb_2_mid_interp(idx_s2);
Vdiff_standard_3 = V_phys_mid_s3 - Vsb_3_mid_interp(idx_s3);

e_standard_2_MAX = max(abs(Vdiff_standard_2)/max(abs(V_phys_mid_s2)));
e_standard_2_AVG = (sum(abs(Vdiff_standard_2))/length(Vdiff_standard_2))/(sum(abs(V_phys_mid_s2))/length(V_phys_mid_s2));
e_standard_3_MAX = max(abs(Vdiff_standard_3)/max(abs(V_phys_mid_s3)));
e_standard_3_AVG = (sum(abs(Vdiff_standard_3))/length(Vdiff_standard_3))/(sum(abs(V_phys_mid_s3))/length(V_phys_mid_s3));

% Error estimate for Taylor basis, only on overlap
V_phys_mid_t2 = V_phys_mid(idx_t2);
V_phys_mid_t3 = V_phys_mid(idx_t3);

Vdiff_taylor_2 = V_phys_mid_t2 - Vtb_2_mid_interp(idx_t2);
Vdiff_taylor_3 = V_phys_mid_t3 - Vtb_3_mid_interp(idx_t3);

e_taylor_2_MAX = max(abs(Vdiff_taylor_2)/max(abs(V_phys_mid_t2)));
e_taylor_2_AVG = (sum(abs(Vdiff_taylor_2))/length(Vdiff_taylor_2))/(sum(abs(V_phys_mid_t2))/length(V_phys_mid_t2));
e_taylor_3_MAX = max(abs(Vdiff_taylor_3)/max(abs(V_phys_mid_t3)));
e_taylor_3_AVG = (sum(abs(Vdiff_taylor_3))/length(Vdiff_taylor_3))/(sum(abs(V_phys_mid_t3))/length(V_phys_mid_t3));


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