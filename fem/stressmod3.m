function [A,S] = stressmod3(omega,upOmega,middleOmega,downOmega,dMid, B1, B2, u)
% Calcule le champs de contrainte a partir d'un champs de deplacement 2D
%
% Parametres:
%   - omega : le maillage support du champs de deplacement
%   - B : le tenseur de Hook dans sa notation de Voigt (matrice)
%   - u : le champs de deplacement 2D
%
% Retourne le champs de contrainte (forme vectorielle)
% [\sigma_11,\sigma_22, \sigma12] en chaque points de Gauss de chaques
% elements

    % Verification des entrees
    assert(isa(omega,'Mesh'),'Objet maillage inconnu');
    assert(isnumeric(B1) && all(size(B1) == [3 3]),'Mauvaise representation du tenseur de Hook');
    assert(isnumeric(u) && numel(u) == 2*omega.nbNodes,'Mauvais champs de deplacement 2D');

    % Integration numerique de Gauss
    [Wg,Xg] = gaussPoints(omega);
    
    Ng = numel(Wg);
    A = zeros(3*omega.nbElems*Ng,1);
    S = zeros(3*omega.nbElems*Ng,1);
    for i=1:omega.nbElems
        ids = omega.elems(i,:); % ids des noeuds
        mapu = bsxfun(@(id,j) (id-1)*2+j,ids(:)',(1:2)'); % index de l'inconnue u
        maps = (i-1)*Ng*3 + (1:Ng*3); % index de l'inconnue sigma
        Xe = omega.nodes(ids,:); % coordonnees de l'element

        M1 = shapesFunctions(omega,Xg,Xe,1);
        test1 = ismember(ids, downOmega.elems, "rows");
        test2 = ismember(ids, middleOmega.elems, "rows");
        test3 = ismember(ids, upOmega.elems, "rows");
        if test1 == 1
            D = kron(eye(Ng),B2);
        elseif (test2 == 1) & (dMid(i)<0)
            D = kron(eye(Ng),B1);
        elseif (test2 == 1) & (dMid(i)>0)
            D = kron(eye(Ng),B1);
        elseif test3 == 1
            D = kron(eye(Ng),B1);
        else
            D = 'error'
        end
        A(maps(:),1) = A(maps(:),1) + D*M1*u(mapu(:)); % Stress (\sigma11, \sigma22, \sigma12)
        S(maps(:),1) = S(maps(:),1) + M1*u(mapu(:)); % Stress (\epsilon11, \epsilon22, 2*\epsion12)
    end
end
