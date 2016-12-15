% example of using local linear regression to find unique eigendirections
% for a swiss roll data set

%% define initial parameters

% number of data points
N = 2000; 

% construct archemedian spiral
a = 1;
theta_vec = linspace(0, 4*pi, 100);
s = 0.5*a*(theta_vec.*sqrt(1+theta_vec.^2)+log(theta_vec+sqrt(1+theta_vec.^2)));

% height
h = 20;

% number of diffusion maps eigenvectors to compute
neigs = 10;

%% generate data

% intialize random number generator
rng(321);

% find angles which correspond to uniform sampling along spiral
theta = interp1(s, theta_vec, rand(N, 1)*max(s));

% data uniformly distributed on swiss roll
global x;
global y;
global z;
z = h*rand(N,1); 
x = a * cos(theta) .* theta;
y = a * sin(theta) .* theta;

% store all data
data = [x y z]; 

%% diffusion maps

% pairwise distance matrix
W = squareform(pdist(data));

% kernel scale
eps = sqrt(5);

% compute embedding coordinates
[V, D] = dmaps(W, eps, neigs);

%% local linear regression

% regression kernel scale
eps_med_scale = 3;

% compute cross-validation error in fitting diffusion maps eigenvectors 
% as a function of previous eigenvectors
res = compute_residuals_DMAPS(V, eps_med_scale);

%% make plots

% plot original data
figure;
plot3(x,y,z,'.')
% view(-20, 80)
xlabel('z_1')
ylabel('z_2')
zlabel('z_3')
axis equal
grid on
title('Original data')

% plot diffusion maps eigenvalues, colored by cross-validation error
figure;
colored_bars(diag(D), res)
set(gca, 'ylim', [0.9 1])
xlabel('k')
ylabel('\mu_k')
axis square
colorbar
title('Diffusion maps eigenvalues, colored by cross-validation error')

% plot data colored by first identified unique eigendirection
figure;
scatter3(x,y,z,50,V(:,2),'.')
% view(-20, 80)
xlabel('z_1')
ylabel('z_2')
zlabel('z_3')
axis equal
grid on
colorbar
title('Data colored by first identified unique eigendirection (k=2)')

% plot data colored by second identified unique eigendirection
figure;
scatter3(x,y,z,50,V(:,5),'.')
% view(-20, 80)
xlabel('z_1')
ylabel('z_2')
zlabel('z_3')
axis equal
grid on
colorbar
title('Data colored by second identified unique eigendirection (k=5)')

% plot data colored by identified repeated eigendirection
figure;
scatter3(x,y,z,50,V(:,3),'.')
% view(-20, 80)
xlabel('z_1')
ylabel('z_2')
zlabel('z_3')
axis equal
grid on
colorbar
title('Data colored by identified repeated eigendirection (k=3)')