function [FX, RES] = local_linear_regression(Y, X, EPS_MED_SCALE)
%LOCAL_LINEAR_REGRESSION Computes the local linear regression fit of a
%function
% [FX, RES] = local_linear_regression(Y, X, EPS_MED_SCALE) computes the
% local linear regression fit of a function Y = f(X)
% 
% Y is an n-long vector of the function output/response at n points
%
% X is an nxk matrix of the k-dimensional function input at n points
%
% EPS_MED_SCALE is the scale to use in the local linear regression kernel 
% the kernel will be a Gaussian with width median(distances)/eps_med_scale
%
% FX is an n-long vector of the local linear approximation to Y at the 
% n points of interest
% 
% RES is an n-long vector of the normalized leave-one-out cross validation 
% error at the n points of interest

% number of data points
[n, k] = size(X);

% local kernel matrix
K = squareform(pdist(X));
eps = median(K(:))/EPS_MED_SCALE;
W = exp(-K.^2 / eps^2);

% compute local fit for each data point
% Storage for alphas (constant terms)
alphas = zeros(n, 1);

% Storage for local linear coefficients.
betas = zeros(n, k);

% Storage for approximation.
FX = zeros(n, 1);

for i=1:n
    % Xx has shape (n, k+1)
    % It concatenates ones with the distance from the focal point i.
    Xi = X-repmat(X(i,:), n, 1);
    Xx = [ones(n, 1)    Xi];

    % Xx2 has shape (k+1, n)
    Xx2 = Xx' * diag(W(i, :));
    
    % Xx2*Xx has shape (k+1, k+1)
    % A has shape (k+1, 1)
    A = (Xx2*Xx)\Xx2 * Y;
    
    % Save fitted coefficients.
    alphas(i) = A(1);
    betas(i, :) = A(2:end);
    
    % Save the approximation.
    FX(i) = alphas(i) + betas(i, :) * X(i, :)';
end


num = 0;
for i=1:n
    num = num + (Y(i) - FX(i))^2;
end

den = 0;
for i=1:n
    den = den + Y(i)^2;
end


RES = sqrt(num / den);
end