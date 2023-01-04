function [A,M] = FEMMatmod(omega,downOmega, B, B2)
% Construit la matrice A elements finis du systeme Au=b, tel que :
%   A = \int_{omega} [\epsilon(N)]^T*B*[\epsilon(N)] dx
%   et optionnellement
%   M = \int_{omega} N^T*eye(2)*N dx
%
%   B est le tenseur de Hook, ecrit sous forme matriciel avec les notations
%   de Voigt.

    % Verification des donnees
    assert(isa(omega,'Mesh'),'Parametre #1 invalide : objet maillage type invalide');
    assert(isnumeric(B) && (all(size(B) == [3 3]) || all(size(B) == [2 2])),'Mauvaise representation du tenseur de Hook');

    A = matrixAssembly(omega,downOmega.elems,middleOmega.elems, 1,B, B2);
    if nargout == 2
        M = matrixAssembly(omega,0,eye(2));
    end
end

function A = matrixAssembly(omega, downOmega,middleOmega, order,B, B2)
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
%keyboard
        [M1,J] = shapesFunctions(omega,Xg,Xe,order); % Evaluation des fonctions de formes

        test = ismember(ids, downOmega, "rows");
        test2 = ismember(ids, middleOmega, 'rows');
        D = kron(diag(Wg.*detJ(J)),B);
        if test(1) == 1
            D = kron(diag(Wg.*detJ(J)),B2);
        end



        if test2(1) == 1 % and the condition for choosing B
            








        A(map(:),map(:)) = A(map(:),map(:))+(M1'*D*M1);
    end
end
