function [x,distG] = signDist(Xe, H, varargin)
%Sign Distance Function
% Parameters = [Xe, H, (N), 'circ, 'r', (N)]
% If number of arguments == 2, it supposes that it is a bimaterial strip
% If 'circ' is used, along with r, as argument, it means there is a
% circular inclusion of radius r, r should be the 4th argument, therefore.
% N is the shape function for finding the interpolated distance,
% with shape function N. Default shape function is N = [0.33,0.33,0.33],
% Separate shape functions can be provided in third or fifth argument
% Note: This takes in only a matrix of three nodes at a time. 

    if nargin==2 || nargin == 3
        x = (Xe(:, 2) - H/2);
    elseif nargin==4 || nargin ==5
        r = varargin{2};
        C = [H/4,H/2];
        x = ((((C(1) - Xe(:,1)).^2)+((C(2) - Xe(:,2)).^2)).^0.5) - r;
    end
    if nargout == 2
        if nargin == 3
            N = varargin{1};
        elseif nargin ==5
            N = varargin{3}; 
        else
            N = [0.3333333; 0.33333333; 1-0.33333333-0.33333333];
        end
        distG = dot(N, x);
    end
end