clear
close all
clc


drag = readmatrix('drag_coefficient.csv');
lift = readmatrix('lift_coefficient.csv');

N = 0:0.125:5;
plot(drag(:,1), drag(:,2), 'r-');
grid on;
title('Drag');
xlabel('Time');
ylabel('Drag coeff value');

figure (2);

plot(lift(:,1), lift(:,2), 'b-');
grid on;
title('Lift');
xlabel('Time');
ylabel('Lift coeff value');