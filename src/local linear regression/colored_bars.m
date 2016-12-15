function colored_bars(BAR_HEIGHT, BAR_COLOR)
%COLORED_BARS makes a colored bar plot
% colored_bars(BAR_HEIGHT, BAR_COLOR) makes a bar plot with variable bar
% height and bar color
%
% BAR_HEIGHT is the height of each bar
%
% BAR_COLOR is the color of each bar (as determined by the corresponding
% colormap); the values of BAR_COLOR are assumed to be between 0 and 1

for i=1:length(BAR_HEIGHT)
    bar(i, BAR_HEIGHT(i), 'facecolor', (1-min(BAR_COLOR(i),1))*[1 1 0]+[0 0 1])
    hold on
end

colormap([linspace(1, 0, 64)' linspace(1, 0, 64)' ones(64, 1)])
caxis([0 1])
set(gca, 'xtick', 1:length(BAR_HEIGHT))