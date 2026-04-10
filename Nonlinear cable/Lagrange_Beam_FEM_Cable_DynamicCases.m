% Non-linear FEM (Modal derivatives - Lagrange) 

clear all
close all

%% Input data

% Material parameters
E = 210e9; % [Pa] Young's modulus
rho = 15.6; % [kg/m] Density of the element
c_damp = 0; % damping coefficient of the element

% Geometric element parameters
d = 50e-3; % [m] diameter of the cable
A = (pi/4)*(d)^2; % [m^2] Cross snnoectional area
I = pi*(d)^4/64; % [m^4] second moment of inertia 
a_L = 1; % [m] length of element
L_tot = 20; % [m] Total length of element (use if easier)


% Load parameters
F1 = -1; % [N] Point load
F2 = -3; % [N/m] Line load
e = 0; % [m] Eccentricity
p = 0;%50; % static line load - set to 0 for cases without and 50 for cases with
P = 4.4237e+03; % Static point load
w1 = (pi/L_tot)^2*sqrt(E*I/rho);
w2 = (2*pi/L_tot)^2*sqrt(E*I/rho);

% Prestressing
N0 = 4.4237e+03;

% Modes and MD's
N_phi = 3;
H = N_phi*(N_phi+1)/2;
MD_phi = 1;

% Switches for different load cases
Loadcase_a = 0;
Loadcase_b = 0;
Loadcase_c = 0;
Loadcase_d = 0;
Loadcase_e = 1;


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

if Loadcase_a == 1 || Loadcase_d == 1
    G = zeros(nel,6);
    for i = 1:nel
      G(i,:) = [E A I rho c_damp 0];
    end
elseif Loadcase_b == 1 || Loadcase_c == 1 || Loadcase_e == 1
    G = zeros(nel,6);
    for i = 1:nel
      G(i,:) = [E A I rho c_damp N0];
    end
end

% understøtninger: U(i)=global dof
U(1) = 1;
U(2) = 2;
U(3) = nno*3-1;

% nodal load (boundary Load): bL(i,:)=[global_dof, magnitude]
%bL = []; %empty if no load
if Loadcase_a == 1 || Loadcase_c == 1
    bL = [nno*3-2 1];
elseif Loadcase_e == 1
    bL = [nno*3-2 P];
else
    bL = [];
end
% element load (domain Load): dL(el,:)=[local_direction, magnitude]
if Loadcase_d == 1
    dL = zeros(nel,2);
    for i = 11:nel
        dL(i,:) = [2 -1];
    end 
elseif Loadcase_e == 1
    dL_static = zeros(nel,2);
    dL_dynamic = zeros(nel,2);
    for i = 11:nel
        dL_static(i,:) = [2 -p];
        dL_dynamic(i,:) = [2 -1];
    end
else
    dL = zeros(nel,2);
end


%% plot geometry
%plotTop(X,T,nel,nno)

%% Constant matrices
dof = 1:3*nno;                         % index to all dofs 
du = U;                                % index to prescribed dofs
df = setdiff(dof,du);                  % index to free dofs

M_mass = M(X,T,G,D,nd,nel);
V = zeros(nd,1);
K0 = Kt(X,T,G,D,nel,nno,V);
K0ff = K0(df,df);                          % Free dofs only 
Mff = M_mass(df,df);                          % Free dofs only


%% Solve linear eigenvalue problem
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

if Loadcase_a == 1 || Loadcase_b == 1 || Loadcase_c == 1
    phi_given = [phi(:,16) phi(:,27) phi(:,34)]; % the 3 first axial modes
    Omega_given = [Omega(16) Omega(27) Omega(34)]; % the corresponding frequencies
else
    phi_given = zeros(size(phi,1),size(phi,2));
    Omega_given = zeros(1,size(phi,2));
end
%% solve complete modal derivative eigenvalue problem

% For the length(df) number of mode-shapes in the system
% and H = length(df)*(length(df)-1)/2 number of modal drivatives in the
% system (when V0=0)

dPhi = zeros(nd,H);
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

%% setting up the modal basis matrix and finding dependent modes

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

%% Program

% Set up time step parameters
dt = 1e-4;
t_CDM = 0:dt:10;

% setup the load vector
if Loadcase_a == 1 || Loadcase_b == 1 || Loadcase_c == 1 || Loadcase_d ==1
    R0 = R(X,T,D,dL,bL,nd,nel); % Load amplitude vector
elseif Loadcase_e == 1
    R0_static = R(X,T,D,dL_static,bL,nd,nel);
    R0_dynamic = R(X,T,D,dL_dynamic,[],nd,nel);
else
    R0 = R(X,T,D,dL,bL,nd,nel);
end
if Loadcase_a == 1 || Loadcase_c == 1
    R_static = P*R0*ones(1,length(t_CDM));
    R_dynamic = zeros(nd,length(t_CDM));
elseif Loadcase_b == 1
    R_static = P*R0*ones(1,length(t_CDM))*0;
elseif Loadcase_d == 1
    R_static = p*R0*ones(1,length(t_CDM));
    t_shift = 3;
    R_dynamic = R0*(F1*sin(w1*(t_CDM-t_shift)) + F2*sin(w2*(t_CDM-t_shift)));
    R_dynamic(:,t_CDM<t_shift) = 0;
elseif Loadcase_e == 1
    R_static = R0_static*ones(1,length(t_CDM));
    t_shift = 3;
    R_dynamic = R0_dynamic*(F1*sin(w1*(t_CDM-t_shift)) + F2*sin(w2*(t_CDM-t_shift)));
    R_dynamic(:,t_CDM<t_shift) = 0;
end

R_combined = R_static + R_dynamic;

if Loadcase_b == 1 || Loadcase_c == 1 || Loadcase_e == 1
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
elseif Loadcase_a == 1 || Loadcase_d == 1
    Q0 = zeros(nd,1);
end

R_load = R_combined - Q0;


% Set up initial state
V0 = zeros(nd,1);
DV0 = zeros(nd,1);
ST0 = zeros(N_phi+H,1);
DST0 = ST0;

% Set up damping
zeta = 0.05;

% Calculate dynamic responses
V_phys = V_physical_pre(X,T,G,D,R_load,U,nno,nel,nd,dt,t_CDM,V0,DV0,zeta,PHI_phys,size(phi,2),Omega);
[Vtb, ST] = V_TaylorBasis_pre(X,T,G,D,R_load,U,nno,nel,nd,N_phi,zeta,dt,t_CDM,ST0,DST0,phi_given,Omega_given);



%% Plot

ux_end_phys = V_phys(nno*3-2,2:end)*1e3;
ux_end_tb = Vtb(nno*3-2,:)*1e3;
uy_mid_phys = V_phys(11*3-1,2:end)*1e3;
uy_mid_tb = Vtb(11*3-1,:)*1e3;

if Loadcase_a == 1 || Loadcase_b == 1 || Loadcase_c == 1
    figure()
    plot(t_CDM,ux_end_phys,'-k')
    hold on
    plot(t_CDM,ux_end_tb,'--k')
    xlabel('time, t [s]')
    ylabel('displacement u_x(L,t)')
    xlim([0 0.3]);
    legend({'Physical u_x(L,t)','Taylor u_x(L,t)'},'Location','best')
elseif Loadcase_d == 1 || Loadcase_e == 1
    figure(1)
    plot(t_CDM,ux_end_phys,':k.')
    hold on
    plot(t_CDM,ux_end_tb,'-k')
    xlabel('time, t [s]')
    ylabel('displacement u_x(L,t)')
    legend({'Physical u_x(L,t)','Taylor u_x(L,t)'},'Location','best')

    figure(2)
    plot(t_CDM,uy_mid_phys,':k.')
    hold on
    plot(t_CDM,uy_mid_tb,'-k')
    xlabel('time, t [s]')
    ylabel('displacement u_y(L/2,t)')
    legend({'Physical u_y(L/2,t)','Taylor u_y(L/2,t)'},'Location','best')
end

if N_phi == 3 && norm(phi_given) == 0 && Loadcase_e == 1

    s1 = ST(1,:);
    s2 = ST(2,:);
    s3 = ST(3,:);
    s4 = ST(4,:);
    s5 = ST(5,:);
    s6 = ST(6,:);
    s7 = ST(7,:);
    s8 = ST(8,:);
    s9 = ST(9,:);
    s4_taylor = s1.^2;
    s5_taylor = s2.*s1;
    s6_taylor = s2.^2;
    s7_taylor = s3.*s1;
    s8_taylor = s3.*s2;
    s9_taylor = s3.^2;

    figure()
    plot(t_CDM,s1,':k.')
    xlabel('time, t [s]')
    ylabel('displacement s')
    legend('s1','Location','best')

    figure()
    plot(t_CDM,s2,':k.')
    xlabel('time, t [s]')
    ylabel('displacement s')
    legend('s2','Location','best')

    figure()
    plot(t_CDM,s3,':k.')
    xlabel('time, t [s]')
    ylabel('displacement s')
    legend('s3','Location','best')

    figure()
    plot(t_CDM,s4,':k.')
    hold on
    plot(t_CDM,s4_taylor,'-k')
    xlabel('time, t [s]')
    ylabel('displacement s')
    legend('s4','s1^2','Location','best')

    figure()
    plot(t_CDM,s5,':k.')
    hold on
    plot(t_CDM,s5_taylor,'-k')
    xlabel('time, t [s]')
    ylabel('displacement s')
    legend('s5','s2*s1','Location','best')

    figure()
    plot(t_CDM,s6,':k.')
    hold on
    plot(t_CDM,s6_taylor,'-k')
    xlabel('time, t [s]')
    ylabel('displacement s')
    legend('s6','s2^2','Location','best')

    figure()
    plot(t_CDM,s7,':k.')
    hold on
    plot(t_CDM,s7_taylor,'-k')
    xlabel('time, t [s]')
    ylabel('displacement s')
    legend('s7','s3*s1','Location','best')

    figure()
    plot(t_CDM,s8,':k.')
    hold on
    plot(t_CDM,s8_taylor,'-k')
    xlabel('time, t [s]')
    ylabel('displacement s')
    legend('s8','s3*s2','Location','best')

    figure()
    plot(t_CDM,s9,':k.')
    hold on
    plot(t_CDM,s9_taylor,'-k')
    xlabel('time, t [s]')
    ylabel('displacement s')
    legend('s9','s3^2','Location','best')

end


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
  s = D*e1 + [N0 0]';
  % make q for local directions
  q = q+w*L*B1'*s;
end

% transform q to global directions
q = A'*q;
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