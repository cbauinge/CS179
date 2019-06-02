clear
clc
close all


M = csvread('build/result.csv');
M = M';

X = linspace(0, 1, size(M, 1));
Y = linspace(0, 1, size(M, 1));

[Xg, Yg] = meshgrid(X, Y);

figure
surf(Xg, Yg, M);
xlabel('x')
ylabel('y')