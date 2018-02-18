function [Cov_X, Cov_XY, Cov_Y] = Cov_comp(X,tau)

%-----------------------------------------------------------------------
% FUNCTION: Cov_comp.m
% PURPOSE:  calculate the covariance matrix of data
% 
% INPUTS:   
%           X: time series data (in form channels x samples)
%           tau: time difference between past and present
% 
% OUTPUT:
%           Cov_X: covariance matrix of X(t-tau) (past)
%           Cov_Y: covariance matrix of X(t) (present)
%           Cov_XY: cross-covariance of X(t-tau) (past state) and X(t) (present state)
%
%-----------------------------------------------------------------------

N = size(X,1);
T = size(X,2);

mstd = 0;
for i=1: N
    mstd = mstd + std(X(i,:));
end
mstd = mstd/N;
X = X/mstd;

t_range1 = 1: 1: T-tau;
t_range2 = 1+tau: 1: T;

X1 = X(:,t_range1);
X1 = bsxfun(@minus, X1, mean(X1,2)); % subtract mean
X2 = X(:,t_range2);
X2 = bsxfun(@minus, X2, mean(X2,2)); % subtract mean

Cov_X = X1*X1'/(T-tau-1); % equal-time covariance matrix at past
Cov_Y = X2*X2'/(T-tau-1); % equal-time covariance matrix at present
Cov_XY = X1*X2'/(T-tau-1); % time-lagged covariance matrix

end