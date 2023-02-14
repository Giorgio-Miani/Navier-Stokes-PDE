clear
close all
clc

ReVStime = readmatrix('');
plot(ReVStime(:, 1), ReVStime(:, 2), ReVStime(:, 3), 'b-');
grid on
title('Reynolds vs time for different meshes');
xlabel('Re');
ylabel('Time');

legend('h = 0.1', 'h = 0.5');

