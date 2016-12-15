function RES = compute_residuals_DMAPS(V, EPS_MED_SCALE)
%COMPUTE_RESIDUALS_DMAPS Computes residuals of local linear fit of
%diffusion maps eigevectors to previous diffusion maps eigenvectors
% RES = compute_residuals_DMAPS(V, EPS_MED_SCALE) computes the local 
% linear regression error for each of the diffusion maps
% eigenvectors as a function of the previous eigenvectors
%
% V are the diffusion maps eigenvectors, stored in columns and ordered by
% the corresponding eigenvalues
% V(:,1) is assumed to be the trivial constant eigenvector
%
% EPS_MED_SCALE is the scale to use in the local linear regression kernel 
% the kernel will be a Gaussian with width median(distances)/eps_med_scale
% (we typically take EPS_MED_SCALE = 3)
% 
% RES are the residuals of each of the fitted functions
% RES(1) is to be ignored (since V(:,1) is the trivial constant 
% eigenvector), and RES(2) will always be 1
% RES(i) is large/close to 1 if V(:,i) parameterizes a new direction in the
% data, and RES(i) is close to 0 if V(:,i) is a harmonic of a previous
% eigenvector

% number of eigenvectors
n = size(V, 2);

% allocate space for residuals
RES = zeros(n,1);

% first residual is 1 by default
RES(2) = 1;

% turn off warnings for singlar matrix
warning('off', 'MATLAB:nearlySingularMatrix');

% compute residuals for local linear fit of eigenvector as a function of
% previous eigenvectors
for i=3:n
    [~, RES(i)] = local_linear_regression(V(:,i), V(:, 2:i-1), EPS_MED_SCALE);
end
