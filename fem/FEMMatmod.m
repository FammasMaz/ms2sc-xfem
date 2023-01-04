function [A,M] = FEMMatmod(omega,middleOmega,downOmega, B, B2, dMid)
% Fork of Chamoin's Code
% Construit la matrice A elements finis du systeme Au=b, tel que :
%   A = \int_{omega} [\epsilon(N)]^T*B*[\epsilon(N)] dx
%   et optionnellement
%   M = \int_{omega} N^T*eye(2)*N dx
% Parameters:  
%   omega, middleOmega, downOmega, can be any matrices defining the mesh of
%   the meshes. (For example, in case of Circular Inclusion, middleOmega,
%   is the mesh around the round material interface, and downOmega is the
%   mesh inside the inclusion. Omega is always the complete mesh.
%   B is the Hook's Tensor (Voigt Notation)

    % Data Verification
    assert(isa(omega,'Mesh'),'Parametre #1 invalide : objet maillage type invalide');
    assert(isnumeric(B) && (all(size(B) == [3 3]) || all(size(B) == [2 2])),'Mauvaise representation du tenseur de Hook');

    A = matrixAssembly(omega,middleOmega, downOmega,1, B, B2, dMid);
    if nargout == 2
        M = matrixAssembly(omega,0,eye(2));
    end
end
function A = matrixAssembly(omega, middleOmega, downOmega,order,B, B2, dMid)
    % Determinant de la matrice jacobienne
    detJ = @(J) J(1:2:end,1).*J(2:2:end,2) - J(2:2:end,1).*J(1:2:end,2);

    % Integration numerique de Gauss
    [Wg,Xg] = gaussPoints(omega,2*(omega.order-order));

    % Creation de la matrice
    A = sparse(2*omega.nbNodes,2*omega.nbNodes);
    for i=1:omega.nbElems
        ids = omega.elems(i,:); % ids des noeuds
        map = bsxfun(@(id,j) (id-1)*2+j,ids(:)',(1:2)'); % index de l'inconnue
        Xe = omega.nodes(ids,:); % coordonnees de l'element
        [M1,J] = shapesFunctions(omega,Xg,Xe,order); % Evaluation des fonctions de formes
        test = ismember(ids, downOmega.elems, "rows"); % Test if element is in downOmega
        test2 = ismember(ids, middleOmega.elems, "rows"); % Test if element is in middleOmega
        
        % Defining B based on location of element
        if test == 1
            D = kron(diag(Wg.*detJ(J)),B2);
        elseif (test2 ==1) & (dMid(i)<0)
            D = kron(diag(Wg.*detJ(J)),B);
        elseif (test2 ==1) & (dMid(i)>0)
            D = kron(diag(Wg.*detJ(J)),B);
        else
            D = kron(diag(Wg.*detJ(J)),B);
        end
        % Constructing Global Matrix of size (nbElements*2 x nbElements*2)
        A(map(:),map(:)) = A(map(:),map(:))+(M1'*D*M1);
    end
end