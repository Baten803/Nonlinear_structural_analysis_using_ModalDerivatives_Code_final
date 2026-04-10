function [V, RR, Rsup, NN, V_storage, R_storage] = NewmarkStatic(X,T,D,G,U,bL,dL,nno,nel,nd)

% set up indexes
dof = 1:nd;                         % index to all dofs 
du = U;                                % index to prescribed dofs
df = setdiff(dof,du);                  % index to free dofs

R0 = R(X,T,D,dL,bL,nd,nel); % Load amplitude vector


%% Newmark method - Static

nstep = 100; %Number of iteration steps
nit = 20; %maximum number of residual updating steps
tol = 1e-5; %displacement vector error tolerance
Vmax = 0.1; %Maximum allowed size of a displacement increment

DRi = R0/nstep; % Base load increment vector (nd x 1)
DR = zeros(nd,1); % Placeholder for initial load step
Ri = zeros(nd,1); % External load vector at current step
V = zeros(nd,1); % Current displacement iterate
dV = zeros(nd,1); % Placeholder for displacement increments
r = zeros(nd,1); % Residual vector

DRpre = zeros(nd,1); % Converged load increment from previous step
DVpre = zeros(nd,1); % Converged displacement increment from previous step

Vpre = zeros(nd,1); % Storerage of converged displacement step
Rpre = zeros(nd,1); % Storage of converged load step

V_storage = zeros(nd,nstep);
R_storage = zeros(nd,nstep);

for i = 1:nstep

  %---------Store converged step from previous step---------------
  Vpre = V;
  Rpre = Ri;

  %---------Start with nominal load increment---------------------
  DR = DRi;

  %---------Initial guess for this steps load vector---------------
  Ri = Rpre + DR;

  %---------Start Newton-Rhapson iteration from previous converged step----
  Vi = Vpre;

  %---------Residual at start of step--------------------------------------
  r = Ri - Q(X,T,G,nel,nno,Vi);
  r(du) = 0;

  %---------Calculate "size" of the residual----------------------
  h = sqrt(r'*r);

  j = 0;
  while h > tol && j < nit
    % Tangent stiffness at current displacement
    Kstar = Kt(X,T,G,nel,nno,Vi);
    
    % Calculate displacement increment
    dV = zeros(nd,1);
    dV(df) = Kstar(df,df)\r(df);
    
    % Update displacement
    Vtrial = Vi + dV;
    
    %------------Length control------------
    DV = Vtrial - Vpre;
    DVnorm = norm(DV);
    if DVnorm > Vmax
      Vtrial = Vpre + (Vmax/DVnorm)*DV;
      DR = (Vmax/DVnorm)*DR;
    end
    
    %-----------Direction control-----------
    Vnorm = DV'*DV;
    Rnorm = DR'*DR;

    dir = (DVpre'*DV)/Vnorm + (DRpre'*DR)/Rnorm;
    if dir < 0
      % change direction of the displacement increment
      Vtrial = Vpre - (Vtrial - Vpre);
      % flip load increment
      DR = -DR;
    end

    % Update the external load with updated load increment
    Ri = Rpre + DR;

    % Accept trial displacement as true displacement
    Vi = Vtrial;

    % Calculate residual based on updated displacement
    r = Ri - Q(X,T,G,nel,nno,Vi);
    r(du) = 0;
    
    % Calculate "size of the updated residual
    h = sqrt(r'*r);

    j = j +1;
  end

  % Accept converged (or last) state for the step
  V = Vi; 

  % Store the accepted displacement and load
  V_storage(:,i) = V;
  R_storage(:,i) = Ri;

  % Store converged increments for direction control
  DVpre = V - Vpre;
  DRpre = Ri - Rpre;

  %Show iteration
  fprintf('step = %d, iterations = %d, ||r|| = %.3e\n', i, j, h);

  if h > tol
    warning('Step %d did not converge. Consider smaller umax or smaller DRshape.', i);
    break
  end

end

%% Reaction forces
% Calculate the reaction forces in the system for every displacement
RR = Q(X,T,G,nel,nno,V); 

% Find the supported dofs
Rsup = RR(du);

%% Normal forces
% Calculate normal forces and moments in each element by use of the 
% generalised strains definition eps = Bv+0.5*I1*v'*G*v
NN = zeros(nel,2);
for nn=1:nel
  no11 = T(nn,1);
  no22 = T(nn,2);
  XX1 = X(no11,:);               % Coordinates of first node in the considered element
  XX2 = X(no22,:);               % Coordinates of second node in the considered element
  Gee = G(nn,:);                    % Material parameters of the given element     
  dee = [no11*3-2 no11*3-1 no11*3 ...
         no22*3-2 no22*3-1 no22*3];
  vee = V(dee);    % Identify global dofs for the given element
  sigma = Sigma(XX1,XX2,Gee,vee);     % Calculate normal forces and moments (the Sigma function is the qbeam function modified very little)
  NN(nn,:) = sigma(1,:);            % Normal force in the end nodes            
end


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


