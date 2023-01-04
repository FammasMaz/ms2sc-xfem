function [A,M] = KaaCirc(inOmega,B1,B2, H,q)
% Construit la matrice A elements finis du systeme Au=b, tel que :
%   A = \int_{omega} [\epsilon(N)]^T*B*[\epsilon(N)] dx
%   et optionnellement
%   M = \int_{omega} N^T*eye(2)*N dx
%
%   B est le tenseur de Hook, ecrit sous forme matriciel avec les notations
%   de Voigt.

    % Verification des donnees
    assert(isa(inOmega,'Mesh'),'Parametre #1 invalide : objet maillage type invalide');
    assert(isnumeric(B1) && (all(size(B1) == [3 3]) || all(size(B1) == [2 2])),'Mauvaise representation du tenseur de Hook');
    oNodes = [1:inOmega.nbNodes];
    [X,Y] = ismember(oNodes,unique(inOmega.elems));
    oNodes(Y(X)) = [];
    A = matrixAssembly(inOmega,1,B1, B2, H,q);
    if nargout == 2
        M = matrixAssembly(inOmega,0,eye(2));
    end
end

function A = matrixAssembly(middleOmega,order,B1, B2, H,q)
    % Determinant de la matrice jacobienne
    detJ = @(J) J(1:2:end,1).*J(2:2:end,2) - J(2:2:end,1).*J(1:2:end,2);

    % Integration numerique de Gauss
    [Wg,Xg] = gaussPoints(middleOmega,2*(middleOmega.order-order));
    a = Xg;% Evaluation des fonctions de formes
    N = [a(:,1) a(:,2) 1-a(:,1)-a(:,2)];
    % Creation de la matrice
    %A = sparse(2*middleOmega.nbNodes,2*middleOmega.nbNodes);
    A = sparse(2*middleOmega.nbNodes,2*(length(unique(middleOmega.elems))));
    nod = unique(middleOmega.elems);
    for i=1:middleOmega.nbElems
        ids = middleOmega.elems(i,:); % ids des noeuds
        idx = [find(nod==ids(1)); find(nod==ids(2)); find(nod==ids(3))]';
        map = bsxfun(@(id,j) (id-1)*2+j,idx(:)',(1:2)'); % index de l'inconnue
        mapr = bsxfun(@(id,j) (id-1)*2+j,ids(:)',(1:2)');% index de l'inconnue
        Xe = middleOmega.nodes(ids,:); % coordonnees de l'element
%keyboard
        [M1,J] = shapesFunctions(middleOmega,Xg,Xe,order);
        DN = [M1(1,1) M1(1, 3) M1(1,5);M1(2,2) M1(2, 4) M1(2,6)];

        if q ==3
            [BI, Bj] = fanB(Xe, H, N, DN,J, B1, B2,'alt');
        else
            [BI, Bj] = fanB(Xe, H, N, DN,J, B1, B2);
        end
        D = kron(diag(Wg.*detJ(J)),Bj);

        A(mapr(:),map(:)) = A(mapr(:),map(:))+(M1'*D*BI);
    end
end
