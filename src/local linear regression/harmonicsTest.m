% Test local linear regression code by explicitly constructing some
% harmonics and non-harmonics.
close all; clear all; clc;


%% Generate test data.
x1 = 0:.01:pi;
x2 = 0:.01:pi;

x1 = x1(randperm(length(x1)));
x2 = x2(randperm(length(x2)));

X = [x1', x2'];

[n, d] = size(X);

frequencies = [0 1 2 3 0 0 0 1 2; ...
               0 0 0 0 1 2 3 2 2]';
[maxk, alsod] = size(frequencies);
           
Y = zeros(n, maxk);
for k=1:maxk
    
    % This works:
    f = frequencies(k, :); f1 = f(1); f2 = f(2);
    % But this doesn't:
    %[f1, f2] = frequencies(k, :)
    % Why? Because MATLAB is terrible.
    
    Y(:, k) = sin(x1 * f1) + sin(x2 * f2);
end


%% Plot the raw data.
% figure();
% subplot(1, 2, 1);
% hold all;
% legends = {};
% for k=1:maxk
%     scatter(x1, Y(:, k), 12, 'filled');
%     f = frequencies(k, :); f1 = f(1); f2 = f(2);
%     legends = [legends, sprintf('f1=%d, f2=%d', f1, f2)];
% end
% legend(legends)
% 
% subplot(1, 2, 2);
% hold all;
% legends = {};
% for k=1:maxk
%     scatter(x2, Y(:, k), 12, 'filled');
%     f = frequencies(k, :); f1 = f(1); f2 = f(2);
%     legends = [legends, sprintf('f1=%d, f2=%d', f1, f2)];
% end
% legend(legends)

% figure();
% for i=1:maxk
%     for j=1:maxk
%         ax = subplot(maxk, maxk, i + (j-1)*maxk);
%         scatter(Y(:, i), Y(:, j), 'k.');
%         set(ax, 'xtick', []);
%         set(ax, 'xticklabel', []);
%         set(ax, 'ytick', []);
%         set(ax, 'yticklabel', []);
%         if i == 1
%             ylabel(sprintf('\\phi_{%d}', j));
%         end
%         
%         if j == 1
%             xlabel(sprintf('\\phi_{%d}', i));
%         end
%         
%         set(gca(), 'XAxisLocation', 'Top');
% %         ax.XAxisLocation = 'origin';
%     end
% end

%% Compute residuals
res = zeros(maxk, 1);
FX = zeros(n, maxk);
res(1) = 0; FX(:, 1) = 0;
res(2) = 1; FX(:, 2) = Y(:, 2);
eps_med_scale = 1;
for k=3:maxk
    [f, r] = local_linear_regression(Y(:, k), Y(:, 2:(k-1)), eps_med_scale);
    res(k) = r;
    FX(:, k) = f;
end
% res = compute_residuals_DMAPS(Y, eps_med_scale)

figure('OuterPosition', [0 0 1600 1000]);
for i=1:maxk
    for j=1:maxk
        ax = subplot(maxk, maxk, i + (j-1)*maxk); hold on;
        scatter(FX(:, i), FX(:, j), 4, 'ro');
        scatter(Y(:, i), Y(:, j), 'k.');

        
        if i == 1
            ylabel(sprintf('f_{%d}', j));
        else
            set(ax, 'ytick', []);
            set(ax, 'yticklabel', []);
        end
        
        
        if j == 1
            xlabel(sprintf('f_{%d}', i));
        else
            set(ax, 'xtick', []);
            set(ax, 'xticklabel', []);
        end
        
        set(gca(), 'XAxisLocation', 'Top');
%         ax.XAxisLocation = 'origin';
    end
end
suptitle('black is truth; red is LLR')
saveas(gcf(), '../../doc/images/harmonicsTest.png');