function [V, D] = dmaps(W, EPS, NCOMPS)
%DMAPS Compute the diffusion maps embedding from the distance matrix W.
% [V, D] = dmaps(W, EPS, NCOMPS) computes the diffusion maps
% embedding
% W is the matrix of pairwise distances
% 
% EPS is the scaling of the diffusion maps kernel
% 
% NCOMPS is the number of embedding coordinates to calculate
% 
% V contains the embedding coordinates (where the i^th column
% contains the i^th emebdding coordinate, and V(:,1) is the trivial 
% constant vector)
%
% D are the eigenvalues correspoinding to the embedding coordinates
    
% compute kernel matrix
W = exp(-(W.^2) ./ (EPS^2));

% compute row sums
d = sum(W);

% normalized symmetric kernel matrix
S = diag(d.^-0.5)*W*diag(d.^-0.5);

% compute eigenvectors
if NCOMPS == max(size(S))
    [V, D] = eig(S);
else
    OPTS.issym = true;
    [V, D] = eigs(S, NCOMPS, 'LM', OPTS);
end

% adjust eigenvectors (since we used symmetrix matrix)
V = diag(d.^-0.5)*V;

% sort eigenvectors and eigenvalues
[~, I] = sort(abs(diag(D)), 'descend');
D = D(I,I);
V = V(:,I);

% normalize eigenvectors
for i=1:NCOMPS
    V(:,i) = V(:,i)./norm(V(:,i));
end

