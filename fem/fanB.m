function [BI, Bj] = fanB(Xe, H, N, DN,J, B1, B2, varargin)
%Local BI, and Bj calculator
%   This function uses is compatible with two distance function
%   formulation, called primary and alternative depending on number of
%   argument.
%   Xe, is nodal coordinates of each element, H is height, r, is radius, N
%   is Shape Function, DN are derivatives of shape functions, J is the
%   Jacobian, B1 and B2 are hook tensors in Voigt Notations, 
%   varargin can contain {'alt'} to solve for alternative formulation
        BI = [];
        [dist, sum] = signDist(Xe,H);
        signD = sign(sum);
        % Use Primary Formulation
        if nargin == 7
            %distN = dot(N,abs(dist));
            distN = abs(dot(N,dist)); % Distance Function
            dSFx = signD*dot(DN(1, :), dist); % Derivative wrt x of Dist fnc
            dSFy = signD*dot(DN(2, :), dist); % Derivative wrt y of Dist fnc 
        %Use Alternative Formulation
        else
            distN = dot(N,abs(dist)) - abs(dot(N, dist)); % Distance Function
            dSFx = dot(DN(1, :), abs(dist))-signD*dot(DN(1, :),dist); % Derivative wrt x of Dist fnc
            dSFy = dot(DN(2, :), abs(dist))-signD*dot(DN(2, :),dist); % Derivative wrt y of Dist fnc
        end
        % Calculate BI for each node
        for j = 1:size(Xe,1)
            locBI = [(DN(1,j)*distN) + (dSFx*N(j)); (DN(2,j)*distN) + (dSFy*N(j))];
            BI(1:2,end+1) = locBI;
        end
        if signD<0 % Check if the element is on either side of discontinuity
            Bj = B2;
        else
            Bj = B1;
        end
        % Reshape BI to be compatible with Voigt Notation
        BI = kron(BI(1:2:end,:),[1 0;0 0;0 1]) + kron(BI(2:2:end,:),[0 0;0 1;1 0]);
end