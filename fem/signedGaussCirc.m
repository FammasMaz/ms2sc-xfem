function [A] = signedGaussCirc(omega, H, r)
%Distance function interpolater for complete meshes
% This function takes in a complete mesh, and for each element,
% interpolates the distance of the Gauss point, using signDist() and maps
% it to the location of the element in the original mesh, contained in A

    [Wg,Xg] = gaussPoints(omega); % Gauss points
    Ng = numel(Wg);
    A = zeros(omega.nbElems*Ng,1);
    N = [Xg(1), Xg(2), 1-Xg(1)-Xg(2)]; % Shape function
    for i=1:omega.nbElems
        maps = (i-1)*Ng*1 + (1:Ng*1); % map for each element
        ids = omega.elems(i,:); % ids des noeuds
        Xe = omega.nodes(ids,:); % coordonnees de l'element
        [x, sum] = signDist(Xe, H, 'circ', r, N); % Caculate distance from Node and Gauss Point
        A(maps(:),1) = A(maps(:),1) + (sum); % Map to correct location of element in omega
    end 
end