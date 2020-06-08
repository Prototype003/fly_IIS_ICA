function [ z ] = fisher_rz(r)
% Converts r-values into z-values
%
% Inputs:
%   r = matrix of r-values
% Outputs:
%   z = matrix of z-values

z = (0.5) * log((1+r)./(1-r));

end

