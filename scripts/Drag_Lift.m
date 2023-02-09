clear
close all
clc

drag = readmatrix('../build/drag_coefficient.csv');
lift = readmatrix('../build/lift_coefficient.csv');

plot(drag(:, 1), drag(:, 2), 'b-');
grid on
title('Drag');
xlabel('Time');
ylabel('Drag');

figure(2);
plot(lift(:,1), lift(:,2), 'b-');
grid on
title('Lift');
xlabel('Time');
ylabel('Lift');