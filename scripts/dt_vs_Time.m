clear
close all
clc

dtVStime = readmatrix('../Results/dt-Tempo.csv');
loglog(1./dtVStime(:, 1), dtVStime(:, 2), 'r--o', 1./dtVStime(:, 1), dtVStime(:, 3), 'b--x', 'LineWidth', 1.5);
hold on
grid on
loglog(1./dtVStime(:, 1), 1./dtVStime(:, 1) .* 500, 'k-.', 'LineWidth', 1.5);
title('dt vs time for different meshes');
xlabel('1/dt');
ylabel('Time');

legend('h = 0.1', 'h = 0.5', '1/dt');