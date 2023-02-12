clear
close all
clc

drag_02 = readmatrix('../Results/Coefficents/h_0.1/Drag/drag_coefficient_Re_2.csv');
drag_05 = readmatrix('../Results/Coefficents/h_0.1/Drag/drag_coefficient_Re_5.csv');
drag_10 = readmatrix('../Results/Coefficents/h_0.1/Drag/drag_coefficient_Re_10.csv');
drag_20= readmatrix('../Results/Coefficents/h_0.1/Drag/drag_coefficient_Re_20.csv');

lift_02 = readmatrix('../Results/Coefficents/h_0.1/Lift/lift_coefficient_Re_2.csv');
lift_05 = readmatrix('../Results/Coefficents/h_0.1/Lift/lift_coefficient_Re_5.csv');
lift_10 = readmatrix('../Results/Coefficents/h_0.1/Lift/lift_coefficient_Re_10.csv');
lift_20= readmatrix('../Results/Coefficents/h_0.1/Lift/lift_coefficient_Re_20.csv');

plot(drag_02(:, 1), drag_02(:, 2), 'k:', 'LineWidth', 1.5);
hold on
grid on
plot(drag_05(:, 1), drag_05(:, 2), 'b-.', 'LineWidth', 1.5);
plot(drag_10(:, 1), drag_10(:, 2),'r--', 'LineWidth', 1.5);
plot(drag_20(:, 1), drag_20(:, 2), 'm-', 'LineWidth', 1.5);
title('Drag coefficents with Re = 2, 5, 10, 20 and h = 0.1');
xlabel('Time');
ylabel('Drag coefficent');
legend('Re = 2', 'Re = 5', 'Re = 10', 'Re = 20');

figure(2)
plot(lift_02(:, 1), lift_02(:, 2), 'k:', 'LineWidth', 1.5);
hold on
grid on
plot(lift_05(:, 1), lift_05(:, 2), 'b-.', 'LineWidth', 1.5);
plot(lift_10(:, 1), lift_10(:, 2),'r--', 'LineWidth', 1.5);
plot(lift_20(:, 1), lift_20(:, 2), 'm-', 'LineWidth', 1.5);
title('Lift coefficents with Re = 2, 5, 10, 20 and h = 0.1');
xlabel('Time');
ylabel('Lift coefficent');
legend('Re = 2', 'Re = 5', 'Re = 10', 'Re = 20');