function [B] = rescale(A, l, u)
%RESCALE Scale range of array elements
% B = rescale(A,l,u) scales the entries of an array to the interval [l,u].
%
% Function exists from MATLAB r2017b onwards
% This function/algorithm is taken from r2017b documentation

inmin = min(A(:));
inmax = max(A(:));

B = l + [(A-inmin)./(inmax-inmin)].*(u-l);

end

