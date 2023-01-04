close all; clear all; clc;
%% Script for Solving XFEM Problem For Circular Discontinuity
% (Q = 2: Primary Distance Formulation, Q = 3: Alternative Distance
% Formulation)
Q = 3; % Question no, in Project Question Paper
%%
addpath('./fem/');
addpath('./meshes/')
addpath('./eet/');
% Parameters of the problem
E1 = 200;
E2 = 100;
nu1 = 0.3;
nu2 = 0.3;
B1 = E1/(1+nu1)/(1-2*nu2)*[1-nu1 nu1 0; nu1 1-nu1 0; 0 0 (1-2*nu1)/2]; % Hook Matrix (With Voigt Notations Voigt)
B2 = E2/(1+nu2)/(1-2*nu2)*[1-nu2 nu2 0; nu2 1-nu2 0; 0 0 (1-2*nu2)/2]; % Hook Matrix (With Voigt Notations Voigt)
Fd = 1000*[0 1]; % Force Definition (x-axis, y-axis)
a = 1;


% Creating Main Meshes
% Mesh, geo2msh, functions taken from Chamoin's Code
omega = Mesh(geo2msh('./meshes/platexfem.geo',1,0.1));

% Creating Secondary Meshes
% Secondary functions taken from Chamoins Code
H = max(omega.nodes(:,2)); % Height of Plate
L = H/2; % Length of the Plate
R = L/4; % Radius of the Inclusion
C = [L/2, H/2]; % Center of the Plate
domega = omega.border; % Borders of the Plate
dtop = domega.restrict(@(x) x(:,2) == H); % Top Edge of the Plate
dbottom = domega.restrict(@(x) x(:,2) == 0); % Bottom Edge of the Plate

% Dividing meshes into three portions: outer, middle, and inner part using
% restrict()
ouOmega = omega.restrict(@(x) ((-x(:,2)+C(2)).^2)+((-x(:,1)+C(1)).^2) >= (R*(1+0.01))^2); % Outer Part
inOmega = omega.restrict(@(x) ((-x(:,2)+C(2)).^2)+((-x(:,1)+C(1)).^2) <= (R*(1-0.01))^2); % Inner Part

% Difference between all elements, and concat of outer and inner part
middleOmega = setdiff(omega.elems, [ouOmega.elems; inOmega.elems], 'rows'); 
% Mesh with the complete number of nodes and middle part elements
middleOmega = Mesh(omega.nodes,middleOmega,omega.type); % Middle Part

dist = signedGaussCirc(middleOmega, H, R); % Dist of middle element wrt Gauss points
Mid = signedGaussCirc(omega, H, R); % Distance of all elements wrt wrt Gauss points
dMid = sign(Mid);% Sign of distances of all elements

% Global stiffness calculation
%K = FEMMatmod(omega,middleOmega,inOmega, B1, B2, dMid);
% Alternative method (Comment one of these)
K = FEMMat(ouOmega,B1) + FEMMat(inOmega,B2) + FEMMat(middleOmega,B1);

% Enrichement Matrices: Kua, Kaa
KUA = Kua_fix(middleOmega, B1, B2, H, Q, 'circ', R);
KAA = Kaa_fix(middleOmega, B1, B2, H, Q, 'circ', R);

% Force definition
F = FEMVec(dtop,Fd); % FEM force vector
Q = zeros(length(KAA),1); % Enrichment force vector
Ftot = [F;Q]; % Augmentation

% Augmented Stiffness
Ktot = [K KUA; KUA' KAA];

% Unknown Displacements (FEM)
u = zeros(size(K,2),1);

% Imposing BC Conditions
% CL() function taken from Chamoin's code
[cl_index,u0] = CL(dbottom, 0, [0 1]); % Imposing 0 displacement in y-axis
cl_index(7) = 1; % Imposing 0 displacement for Left Bottom Edge
u0 = [u0;0]; % Augmenting u0 for extra BC
u(cl_index) = u0; % Predefining displacements in global u
utot = [u;Q]; % Augmenting u with initial enriched displacements 

cl_index = [cl_index; logical(Q)]; % Augmenting cl_index with enriched nodes indices 

% Solve the augmented system
utot(~cl_index) = Ktot(~cl_index,~cl_index)\(Ftot(~cl_index) - Ktot(~cl_index,cl_index)*utot(cl_index));

% Extract the XFEM Displacements
u = utot(1:size(K,1));

% Calculation of stress and strain
[sigma, strains] = stressmod3(omega,ouOmega, middleOmega,inOmega, dMid, B1,B2,u);
%%
e1 = strain(inOmega,u);
e2 = strain(ouOmega,u);
s1 = stress(inOmega, B2, u);
s2 = stress(ouOmega, B1, u);

% Plot Stress, Displacement, and Strain Fields
figure('Name','Solution');
  subplot(1,3,1)
    plotElemField(deform(omega,u,1./max(abs(u))),sigma);
    xlabel('x');
    ylabel('y');
    hc=colorbar;
    title(hc,'Pascals');
    title('Stresses');
  subplot(1,3,2);
    plot(omega.border);
    hold on
    quiver(omega.nodes(:,1),omega.nodes(:,2),u(1:2:end),u(2:2:end));
    xlabel('x');
    ylabel('y');
    title('Displacement');
  subplot(1,3,3);
    plotElemField(deform(omega,u,1./max(abs(u))),strains);
    xlabel('x');
    ylabel('y');
    hc = colorbar;
    title(hc,'Unitless');
    title('Strain');

 figure('Name','Strains');
  subplot(1,2,1)
    plotElemField(deform(inOmega,u,1./max(abs(u))),e1);
    plotElemField(deform(ouOmega,u,1./max(abs(u))),e2);
    xlabel('x');
    ylabel('y');
    hc=colorbar;
    title(hc,'Unitless');
    title('Strains');
  subplot(1,2,2);
    plotElemField(deform(inOmega,u,1./max(abs(u))),s1);
    plotElemField(deform(ouOmega,u,1./max(abs(u))),s2);
    xlabel('x');
    ylabel('y');
    hc = colorbar;
    title(hc,'Pascals');
    title('Stresses');
