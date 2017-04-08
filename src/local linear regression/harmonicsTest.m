% Test local linear regression code by explicitly constructing some
% harmonics and non-harmonics.
%close all; clear all; clc;
if 0
    methodString = 'originalResidualCode';
    titleString = 'original code';
    method = @(Ythis, Yprev, EPS) local_linear_regression_old(Ythis, Yprev, EPS);
else
    methodString = 'rewrittenResidualCode';
    titleString = 'rewritten code';
    method = @(Ythis, Yprev, EPS) local_linear_regression(Ythis, Yprev, EPS);
end
rng(4);


%% Generate test data.
noiseLevel = 0.05;
x1 = 0:.01:pi;
% x1 = x1 + randn(size(x1))*noiseLevel;
x2 = 0:.01:pi;
% x2 = x2 + randn(size(x1))*noiseLevel; 

x1 = x1(randperm(length(x1)));
x2 = x2(randperm(length(x2)));

X = [x1', x2'];

[n, d] = size(X);

frequencies = [0 1 2 3 0 0 0 1 2 4 5 6 7 8 9 4 5 6 7 8 9; ...
               0 0 0 0 1 2 3 2 2 0 0 0 0 0 0 0 0 0 0 0 0]';
[maxk, alsod] = size(frequencies);
           
Y = zeros(n, maxk);
for k=1:maxk
    
    % This works:
    f = frequencies(k, :); f1 = f(1); f2 = f(2);
    % But this doesn't:
    %[f1, f2] = frequencies(k, :)
    % Why? Because MATLAB is terrible.
    
    Y(:, k) = cos(x1 * f1) + cos(x2 * f2) + randn(size(x1)) * noiseLevel;
end


% %% Plot the raw data.
% figure('OuterPosition', [0 0 800 500]);
% subplot(1, 2, 1);
% hold all;
% legends = {};
% for k=1:5
%     scatter(x1, Y(:, k)', 12, 'filled');
%     f = frequencies(k, :); f1 = f(1); f2 = f(2);
%     legends = [legends, sprintf('m=%d, n=%d', f1, f2)];
% end
% xlabel('x')
% ylabel('f_k(x, y)')
% legend(legends, 'Location', 'south')
% 
% subplot(1, 2, 2);
% hold all;
% legends = {};
% for k=1:maxk
%     scatter(x2, Y(:, k)', 12, 'filled');
%     f = frequencies(k, :); f1 = f(1); f2 = f(2);
%     legends = [legends, sprintf('m=%d, n=%d', f1, f2)];
% end
% xlabel('y')
% ylabel('f_k(x, y)')
% legend(legends, 'Location', 'south')
% 
% saveas(gcf(), '../../doc/images/harmonicsTest-data.png');


%% Compute residuals
res = zeros(maxk, 1);
FX = zeros(n, maxk);
res(1) = 0; FX(:, 1) = 2;
res(2) = 1; FX(:, 2) = Y(:, 2);
eps_med_scale = 10;
for k=3:maxk
    [f, r] = method(Y(:, k), Y(:, 2:(k-1)), eps_med_scale);
    res(k) = r;
    FX(:, k) = f;
end
% res = compute_residuals_DMAPS(Y, eps_med_scale)

% figure('OuterPosition', [0 0 1600 1000]);
% for i=1:maxk
%     for j=1:maxk
%         ax = subplot(maxk, maxk, i + (j-1)*maxk); hold on;
%         scatter(FX(:, i), FX(:, j), 4, 'ro');
%         scatter(Y(:, i), Y(:, j), 'k.');
%         xlim([-2, 2]);
%         ylim([-2, 2]);
%         
%         if i == 1
%             ylabel(sprintf('f_{%d}', j));
%         else
%             set(ax, 'ytick', []);
%             set(ax, 'yticklabel', []);
%         end
% 
%         
%         if j == 1
%             xlabel(sprintf('f_{%d}', i));
%         else
%             set(ax, 'xtick', []);
%             set(ax, 'xticklabel', []);
%         end
%         
%         set(gca(), 'XAxisLocation', 'Top');
% %         ax.XAxisLocation = 'origin';
%     end
% end
% suptitle(sprintf('black is truth; red is LLR; %s', titleString));
% saveas(gcf(), sprintf('../../doc/images/harmonicsTest-%s.png', methodString));

figure();
colored_bars(1 ./ (1:maxk), res);
xlabel('k');
ylabel('\lambda_k');
suptitle(sprintf('fake spectrum colored by r_k; %s', titleString));
saveas(gcf(), sprintf('../../doc/images/harmonicsTest-bars-%s.png', methodString));
