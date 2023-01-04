close all; clear all; clc;
%% Script d'exemple de resolution du probleme elements finis
addpath('./fem/');
addpath('./eet/');
%keyboard
% Parametres du probleme
E1 = 200;
E2 = 100;
nu1 = 0.3;
nu2 = 0.3;
B1 = E1/(1+nu1)/(1-2*nu2)*[1-nu1 nu1 0; nu1 1-nu1 0; 0 0 (1-2*nu1)/2]; % Matrix de Hook (Avec les notations de Voigt)
B2 = E2/(1+nu2)/(1-2*nu2)*[1-nu2 nu2 0; nu2 1-nu2 0; 0 0 (1-2*nu2)/2]; % Matrix de Hook (Avec les notations de Voigt)
Fd = 1000*[0 1];
a = 1;

% Charge le maillage
omega = Mesh(geo2msh('./meshes/platecircle.geo',1,0.5));
%keyboard
% Extraction de maillage secondaires
H = max(omega.nodes(:,2));
L = H/2; % Length of the Plate
R = L/4; % Radius of the Inclusion
C = [L/2, H/2]; % Center of the Plate
domega = omega.border;
dtop = domega.restrict(@(x) x(:,2) == H);
dbottom = domega.restrict(@(x) x(:,2) == 0);

ouOmega = omega.restrict(@(x) ((-x(:,2)+C(2)).^2)+((-x(:,1)+C(1)).^2) >= (R-0.01)^2); % Outer Part% Difference between all elements, and concat of outer and inner part
inOmega = omega.restrict(@(x) ((-x(:,2)+C(2)).^2)+((-x(:,1)+C(1)).^2) <= (R+0.01)^2); % Outer Part% Difference between all elements, and concat of outer and inner part

%%

% Ecriture du systeme
K = FEMMat(ouOmega, B1) + FEMMat(inOmega, B2); 
F = FEMVec(dtop,Fd);

% Initialisation de l'inconnue
u = zeros(size(K,2),1);

% Imposition des conditions aux limites
[cl_index,u0] = CL(dbottom, 0, [0 1]);
cl_index(7) = 1;
u0 = [u0;0];
u(cl_index) = u0;



   
% Resolution
u(~cl_index) = K(~cl_index,~cl_index)\(F(~cl_index) - K(~cl_index,cl_index)*u(cl_index));
[sigma, strain] = stressmod(omega,inOmega, B1,B2,u);

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
        plotElemField(deform(omega,u,1./max(abs(u))),strain);
        xlabel('x');
        ylabel('y');
        hc = colorbar;
        title(hc,'Unitless');
        title('Strain');
