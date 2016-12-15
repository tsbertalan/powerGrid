% example of using local linear regression to find unique eigendirections
% for a swiss roll data set
clear all;
close all;

%% Load datafile.
load('../../data/generatorLoads.mat');
[N, ndims] = size(data);

% Normalize the data.
for j=1:ndims
    data(:, j) = (data(:, j) - mean(data(:, j))) / std(data(:, j));
end

% Arbitrarily say that the last few hours' load values are our 'x', 'y',
% and 'z' features for plotting purposes.
fromend = 12;
xinds = [1, size(data, 2)/2, size(data, 2)];
x = data(:, xinds(1));
y = data(:, xinds(2));
z = data(:, xinds(3));

figure(); plot(data', 'k-');
xlim([1, size(data, 2)]);
xlabel('time [h]');
ylabel('generator load');
title('raw data');

%% diffusion maps

% pairwise distance matrix
dz = squareform(pdist(data));

% Try many epsilons.
% compute kernel matrix
epsVals = logspace(-3, 3, 64);
ksumVals = [];
for EPS=epsVals
    W = exp(-(dz.^2) ./ (EPS^2));

    % compute row sums
    d = sum(W);

    % normalized symmetric kernel matrix
    normPow = 1.0;
    S = diag(d.^-normPow)*W*diag(d.^-normPow);
    
    ksumVals = [ksumVals; sum(sum(W))];
end

[~, imax] = max(diff(ksumVals));
epsMax = epsVals(imax);
eps = 4.0;


figure();
semilogx(epsVals, ksumVals, 'k--');
line([eps eps], ylim(), 'Color', 'black');
line([epsMax epsMax], ylim(), 'Color', 'red');
line([median(dz(:)), median(dz(:))], ylim(), 'Color', 'blue');
legend('sum', 'Chosen \epsilon', 'maximum slope \epsilon', 'median distance',...
    'Location', 'west');
xlabel('\epsilon');
ylabel('\Sigma_{ij} exp(-||z_i - z_j||^2/\epsilon^2 )');


%% Actually do Dmaps.
% We certainly don't want to compute (too many) more output features than we have
% input features.
neigs = ndims;

% compute embedding coordinates
[V, D] = dmaps(dz, eps, neigs);

%% local linear regression

% regression kernel scale
eps_med_scale = .5;

% compute cross-validation error in fitting diffusion maps eigenvectors 
% as a function of previous eigenvectors
res = compute_residuals_DMAPS(V, eps_med_scale);
res'

%% Decide on first two best coordinates.
[~, coordOrder] = sort(res, 'descend');
k1 = coordOrder(1);
k2 = coordOrder(2);
for k=2:numel(res)
    if res(k) < .5
        harm1 = k;
        break;
    end
end

%% plot diffusion maps eigenvalues, colored by cross-validation error
figure();
colored_bars(sqrt(diag(D)), res)
set(gca, 'ylim', [0.0 1])
xlabel('k')
ylabel('$\sqrt{\mu_k}$', 'Interpreter', 'latex')

axis square
colorbar
title('Diffusion maps eigenvalues, colored by cross-validation error')

%% 3d plots
% plot original data
figure('Position', [10, 10, 1000, 1000]);
subplot(2, 2, 1);
plot3(x,y,z,'.')
% view(-20, 80)
xlabel(sprintf('load(t=%d)', xinds(1))); ylabel(sprintf('load(t=%d)', xinds(2))); zlabel(sprintf('load(t=%d)', xinds(3)));
axis equal
grid on
title('Original data')

% plot data colored by first identified unique eigendirection
subplot(2, 2, 2);
scatter3(x,y,z,50,V(:, k1),'.')
% view(-20, 80)
xlabel(sprintf('load(t=%d)', xinds(1))); ylabel(sprintf('load(t=%d)', xinds(2))); zlabel(sprintf('load(t=%d)', xinds(3)));
axis equal
grid on
colorbar
title(sprintf('first identified unique\neigendirection (k=%d)', k1))

% plot data colored by second identified unique eigendirection
subplot(2, 2, 3);
scatter3(x,y,z,50,V(:, k2),'.')
% view(-20, 80)
xlabel(sprintf('load(t=%d)', xinds(1))); ylabel(sprintf('load(t=%d)', xinds(2))); zlabel(sprintf('load(t=%d)', xinds(3)));
axis equal
grid on
colorbar
title(sprintf('second identified unique\neigendirection (k=%d)', k2))

% plot data colored by identified repeated eigendirection
subplot(2, 2, 4);
scatter3(x,y,z,50,V(:, harm1),'.')
% view(-20, 80)
xlabel(sprintf('load(t=%d)', xinds(1))); ylabel(sprintf('load(t=%d)', xinds(2))); zlabel(sprintf('load(t=%d)', xinds(3)));
axis equal
grid on
colorbar

title(sprintf('first identified repeated\neigendirection (k=%d)', harm1))