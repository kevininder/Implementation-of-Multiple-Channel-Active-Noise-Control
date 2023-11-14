% Performance evaluation of 1x2x2.

close all; clear; clc;

% read the data gathered in the main file
tbl = readtable('1x2x2.xlsx');

% plot color scale bubble chart
tiledlayout(1, 1)
nexttile
bubblechart(tbl, 'Decreased_decibels', 'Theta_between_error_mics', 'Area', 1, 'MarkerEdgeColor', 'k');
colorbar('Ticks', [9,18,27,36,45], ...
         'TickLabels', {'Ineffective', 'Low', 'Average', 'Good', 'Excellent'}, 'FontSize', 12)

title('1x2x2', 'FontSize', 12)
blgd = bubblelegend('Area(m^{2})', 'FontSize', 12);
lgd = legend('Noise: 70(dB)\newlineAverage required time: 0.0167(s)', 'FontSize', 12);
blgd.Layout.Tile = 'east';
lgd.Layout.Tile = 'east';