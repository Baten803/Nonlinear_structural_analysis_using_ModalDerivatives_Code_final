% Non-linear FEM (Lagrange element) 

clear all
close all

%% Input data

% Material parameters
E_deck = 210e9; % [Pa] Young's modulus
E_hanger = 195e9;
E_cable = 195e9;
rho = 7850; % [kg/m^3] Density of the element
c_damp = 0; % damping coefficient of the element

% Geometric element parameters
A_deck = 1.24; % [m^2] Cross sectional area
A_hanger = 2*0.0162;
A_cable = 2*0.389;
I_deck = 3.94; % [m^4] second moment of inertia
I_hanger = 2*2.09e-5;
I_cable = 2*0.012;
a_L = 1; % [m] length of element
L_tot = 78; % [m] Total length of horizontal structure (use if easier)
h1 = 4; % [m] height of cable at midpoint
h2 = 5;
H = 6.5; % [m] Height of cable supports
Lc1 = sqrt((L_tot/4)^2 + (h2-h1)^2); % Length of inclined cables
Lc2 = sqrt((L_tot/4)^2 + (H-h2)^2);

% Load parameters
p_dead = rho*A_deck*9.81; %[N/m] self-weight

% Prestressing
N0_hangers = 1.23*p_dead*L_tot/4;
N0_cable_inner = N0_hangers*(1/2)*(Lc1/(h2-h1));
N0_cable_outer = N0_hangers*(3/2)*(Lc2/(H-h2));

% Mode parameters
N_phi = 9;
HH = N_phi*(N_phi+1)/2;
MD_phi = 1;

%% Matrices for the element 

% knudekoordinater: X(knude,:)=[x,y] 
Numb = 20; % Number of elements wanted in a straight line
NumbC = 1; % elements of a cable

X = zeros(Numb+1+4*NumbC+1 ,2);

% Horizontal beam nodes
for i = 1:Numb+1
  X(i,:) = [(i-1)*L_tot/Numb 0]; 
end

% Inclined outer most left cable nodes
idx_incleftnodes1 = zeros(NumbC+1,1);
for i = 1:NumbC+1
  idx = i + (Numb+1);
  X(idx,:) = [(i-1)*L_tot/4 H-4*((H-h2)/L_tot)*(i-1)*L_tot/4];
  idx_incleftnodes1(i) = idx;
end

% Inclined inner left cable nodes
idx_incleftnodes2 = zeros(NumbC,1);
idx = Numb+NumbC+3;
X(idx,:) = [2*L_tot/4 h1];
idx_incleftnodes2(1) = idx;

% Inclined inner right cable nodes
idx_incrightnodes1 = zeros(NumbC,1);
idx = Numb+2*NumbC+3;
X(idx,:) = [3*L_tot/4 h2];
idx_incrightnodes1(1) = idx;

% Inclined outer most right cable nodes
idx_incrightnodes2 = zeros(NumbC,1);
idx = Numb+3*NumbC+3;
X(idx,:) = [L_tot H];
idx_incrightnodes2(1) = idx;

nno = size(X,1);    % number of nodes, antal knuder i systemet

% elementer: T(el,:)=[startknude slutknude]

T = zeros(Numb+4*NumbC+3,2);

% Horizontal beam elements
Tflat = zeros(Numb,1);
for i = 1:Numb
  T(i,:) = [i i+1];
  Tflat(i) = i;
end

% Inclined cable elements
Tinc = zeros(2*NumbC,1);
for i = 1:4*NumbC
  T(i+Numb,:) = [i+(Numb+1) i+1+(Numb+1)];
  Tinc(i) = i+Numb; 
end
Tincleft1 = Tinc(1);
Tincleft2 = Tinc(2);
Tincright1 = Tinc(3);
Tincright2 = Tinc(4);

% Vertical cable elements - midpoint
T(Numb+4*NumbC+1,:) = [Numb/2+1 Numb+1+3];
Tv = Numb+4*NumbC+1;

% Vertical cable element - Left
T(Numb+4*NumbC+2,:) = [Numb/4+1 Numb+1+2];
Tl = Numb+4*NumbC+2;

% Vertical element - right
T(Numb+4*NumbC+3,:) = [3*Numb/4+1 Numb+1+4];
Th = Numb+4*NumbC+3;

nel = size(T,1);    % number of elements, antal elementer i systemet

% frihedsgrader, globale dof: D(el,:)=[V1 V2 V3 V4 V5 V6]

D = zeros(nel,6);

for i = 1:nel
  D(i,:) = [3*T(i,1)-2 3*T(i,1)-1 3*T(i,1) 3*T(i,2)-2 3*T(i,2)-1 3*T(i,2)];
end
% % Hinge at the point where all 3 cables meet
% D(31,3) = nno*3+1;
% D(50,6) = nno*3+2;
 
nd=max(max(D));     % number of dofs, antal frihedsgrader

% materialer: G(el,:)=[E-modul, tværsnitsareal, inertimoment]

G = zeros(nel,6);
for i = 1:nel
  if ismember(i,Tflat) 
    G(i,:) = [E_deck A_deck I_deck rho c_damp 0];
  elseif ismember(i,Tincleft1)
    G(i,:) = [E_cable A_cable I_cable rho c_damp N0_cable_outer];
  elseif ismember(i,Tincleft2)
    G(i,:) = [E_cable A_cable I_cable rho c_damp N0_cable_inner];
  elseif ismember(i,Tincright1)
    G(i,:) = [E_cable A_cable I_cable rho c_damp N0_cable_inner];
  elseif ismember(i,Tincright2)
    G(i,:) = [E_cable A_cable I_cable rho c_damp N0_cable_outer];
  elseif ismember(i,Tl)
    G(i,:) = [E_hanger A_hanger I_hanger rho c_damp N0_hangers];
  elseif ismember(i,Tv)
    G(i,:) = [E_hanger A_hanger I_hanger rho c_damp N0_hangers];
  elseif ismember(i,Th)
    G(i,:) = [E_hanger A_hanger I_hanger rho c_damp N0_hangers];
  end
end


% understøtninger: U(i)=global dof
U(1) = 1;
U(2) = 2;
U(3) = (Numb+1)*3-2;
U(4) = (Numb+1)*3-1;
U(5) = (Numb+2)*3-2;
U(6) = (Numb+2)*3-1;
U(7) = (Numb+1+4*NumbC+1)*3-2;
U(8) = (Numb+1+4*NumbC+1)*3-1;

% nodal load (boundary Load): bL(i,:)=[global_dof, magnitude]
bL = []; %empty if no load
%bL = [nno*3-1 -P];

% element load (domain Load): dL(el,:)=[local_direction, magnitude]
dL = zeros(nel,2);
for i = 1:Numb
    dL(i,:) = [2 -p_dead];
end

%% plot geometry
plotTop(X,T,nel,nno)

%% Program

%% Solving the linearized eigenvalue problem

dof = 1:3*nno;                         % index to all dofs 
du = U;                                % index to prescribed dofs
df = setdiff(dof,du);                  % index to free dofs

V = zeros(nd,1);
K0 = Kt(X,T,G,D,nel,nno,V);
M_mass = M(X,T,G,D,nd,nel);

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

m_phi = phi'*M_mass*phi; % making the modal matrices based on the new mode-shapes

%% Finding the modal derivatives

dPhi = zeros(nd,HH);
for i = 1:size(phi,2)
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

% Initialize the PHI-matrix containing the modes and the modal derivatives
PHI_phys = [phi dPhi];
PHi_phys = [phi(df,:) dPhi(df,:)];

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

%% Calculating the dynamic response

%Set up the load vector
R_static = R(X,T,G,D,dL,bL,nd,nel);
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
R_prestress = R_static - Q0;

% Calculate the displacements, reaction forces and normal force in beam
% elements
V0 = zeros(nd,1);
DV0 = zeros(nd,1);
ST0 = zeros(N_phi+HH,1);
DST0 = zeros(N_phi+HH,1);
phi_given = 0;
Omega_given = 0;
zeta = 0.05;
dt_phys = 1e-5;
dt_tb = 1e-4;
t_CDM_phys = 0:dt_phys:30;
t_CDM_tb = 0:dt_tb:30;
R_dynamic_phys = R_prestress*ones(1,length(t_CDM_phys));
R_dynamic_tb = R_prestress*ones(1,length(t_CDM_tb));
%V_phys = V_physical_pre(X,T,G,D,R_dynamic_phys,U,nno,nel,nd,dt_phys,t_CDM_phys,V0,DV0,zeta,PHI_phys,size(phi,2),Omega);
[V_tb, ST_CDM] = V_TaylorBasis_pre(X,T,G,D,R_dynamic_tb,U,nno,nel,nd,N_phi,zeta,dt_tb,t_CDM_tb,ST0,DST0,phi_given,Omega_given);

%V_phys_mid = V_phys(11*3-1,2:end);
V_tb_mid = V_tb(11*3-1,:);

%% Plot
figure()
%plot(t_CDM_phys,V_phys_mid*1e3,':k.')
% hold on 
plot(t_CDM_tb,V_tb_mid*1e3,'-k')
xlabel('Time, t [s]')
ylabel('Vertical displacement u_y(L/2,t)')
legend('physical u_y(L/2,t)','Taylor u_y(L/2,t)','Location','best')










%% functions

function Q1 = Q(X,T,G,D,nel,nno,V)
% Setup system force vector Q
  Q1 = zeros(nno*3,1);                   % initialization of Q1
  for el = 1:nel  
    no1 = T(el,1);  no2 = T(el,2);       % start node/end node
    X1 = X(no1,:);  X2 = X(no2,:);       % coordinates to start/end node
    de = D(el,:);                        % dofs for element
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

function sigma = Sigma(X1,X2,Ge,ve)
% Setup element forces N and M in the end-points of the element

% nip = 3; % Can be used up to polynomial order = 2*nip-1 = 5
% zp = sqrt(3/5);
% zip = [1-zp 1 1+zp]/2;
% wip = [5/9 8/9 5/9]/2;

EA = Ge(1)*Ge(2);
EI = Ge(1)*Ge(3);
D1 = [EA 0; 0 EI];                      % material stiffness matrix
I1 = [1;0];
[A, L] = Abeam(X1,X2);                 % transformation matrix
vel = A*ve;                            % local displacement vector
%q = zeros(6,1);
sigma = zeros(2,2);     

for ii=1:2
  z = (ii-1);
% w = wip(i);
  B = Bbeam(L,z);
  G1 = Gbeam(L,z);
% B1 = B+I1*(vel'*G);
  % strains
  e1 = B*vel+1/2*I1*(vel'*G1*vel);

  % stresses 
  sigma(:,ii) = D1*e1;
% % make q for local directions
% q = q+w*L*B1'*s;
end
% transform q to global directions
%q = A'*q;
end

function K = K(X,T,G,D,nel,nno,V)
% Setup system stiffnessmatrix
K = zeros(nno*3,nno*3);                % initialization of K
for el = 1:nel
  no1 = T(el,1);  no2 = T(el,2);       % start node/end node
  X1 = X(no1,:);  X2 = X(no2,:);       % coordinates to start/end node
  de = D(el,:);                        % dofs for element
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

function Kt = Kt(X,T,G,D,nel,nno,V)
% Setup system stiffnessmatrix
Kt = zeros(nno*3,nno*3);                % initialization of K
for el = 1:nel
  no1 = T(el,1);  no2 = T(el,2);       % start node/end node
  X1 = X(no1,:);  X2 = X(no2,:);       % coordinates to start/end node
  de = D(el,:);                        % dofs for element
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
  s = EA*e1(1) + N0;
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
function R = R(X,T,G,D,dL,bL,nd,nel)
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
R = R;
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
