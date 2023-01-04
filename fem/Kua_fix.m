function [A,M] = Kua_fix(middleOmega,B1,B2, H,q, varargin)
% Construit la matrice A elements finis du systeme Au=b, tel que :
%   A = \int_{omega} [\epsilon(N)]^T*B*[\epsilon(N)] dx
%   et optionnellement
%   M = \int_{omega} N^T*eye(2)*N dx
%
%   B1 and B2 are the Hook tensors, in matrix form, based on Voigt Notation
%   H, is the height, q refers to the distance function formulation (2 =
%   primary, 3 = alternative)
%   varargin can contain {'circ', r} pointing to solving circular inclusion
%   problem

    % Data assertion
    assert(isa(middleOmega,'Mesh'),'Parametre #1 invalide : objet maillage type invalide');
    assert(isnumeric(B1) && (all(size(B1) == [3 3]) || all(size(B1) == [2 2])),'Mauvaise representation du tenseur de Hook');
    
    % Horizontal Problem
    if nargin == 5 
        A = matrixAssembly(middleOmega,1,B1, B2, H,q);
    % Circular Inclusion Problem
    else 
        R = varargin{2}; % Radius
        A = matrixAssembly(middleOmega,1,B1, B2, H,q,R); 
    end
    if nargout == 2
        M = matrixAssembly(middleOmega,0,eye(2));
    end
end

function A = matrixAssembly(middleOmega,order,B1, B2, H,q, varargin)
    
    % Determinant de la matrice jacobienne
    detJ = @(J) J(1:2:end,1).*J(2:2:end,2) - J(2:2:end,1).*J(1:2:end,2);

    % Integration numerique de Gauss
    [Wg,Xg] = gaussPoints(middleOmega,2*(middleOmega.order-order));
    a = Xg;% Evaluation des fonctions de formes
    N = [a(:,1) a(:,2) 1-a(:,1)-a(:,2)];
    % Creation de la matrice
    %A = sparse(2*middleOmega.nbNodes,2*middleOmega.nbNodes);
    A = sparse(2*middleOmega.nbNodes,2*(length(unique(middleOmega.elems))));
    nod = unique(middleOmega.elems); % Enriched Nodes
    for i=1:middleOmega.nbElems
        ids = middleOmega.elems(i,:); % ids des noeuds
        idx = [find(nod==ids(1)); find(nod==ids(2)); find(nod==ids(3))]';
        map = bsxfun(@(id,j) (id-1)*2+j,idx(:)',(1:2)'); % Map for enriched elements
        mapr = bsxfun(@(id,j) (id-1)*2+j,ids(:)',(1:2)');% Map for elements in omega
        Xe = middleOmega.nodes(ids,:); % Coordinates of Elements
%keyboard
        [M1,J] = shapesFunctions(middleOmega,Xg,Xe,order); % Derivative of Shp Fnc and Jacobian
        DN = [M1(1,1) M1(1, 3) M1(1,5);M1(2,2) M1(2, 4) M1(2,6)];
        % Solve Horizontal Problem
        if nargin == 6
            if q ==3
                [BI, Bj] = fanB(Xe, H, N, DN,J, B1, B2,'alt');
            else
                [BI, Bj] = fanB(Xe, H, N, DN,J, B1, B2);
            end
        % Solve Circular Inclusion Problem
        else
            r = varargin{1};
            if q ==3
                [BI, Bj] = fanCirc(Xe, H, r, N, DN,J, B1,B2,'alt');
            else
                [BI, Bj] = fanCirc(Xe, H, r, N, DN,J, B1,B2);
            end
        end
        D = kron(diag(Wg.*detJ(J)),Bj); % Transform D from reference to global

        A(mapr(:),map(:)) = A(mapr(:),map(:))+(M1'*D*BI); % Map Stiffness to right Position
    end
end
