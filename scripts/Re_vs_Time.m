clear
close all
clc

ReVStime = readmatrix('../Results/Reynolds-Tempo.xlsx');
plot(ReVStime(:, 1), ReVStime(:, 2), 'r--o', ReVStime(:, 1), ReVStime(:, 3), 'b--x', 'LineWidth', 1.5);
grid on
title('Reynolds vs time for different meshes');
xlabel('Re');
ylabel('Time');

legend('h = 0.1', 'h = 0.5');

