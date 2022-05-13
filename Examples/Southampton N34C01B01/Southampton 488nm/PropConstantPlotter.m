% Ralf Mouthaan
% University of Cambridge
% December 2021
%
% Script to plot propagation constant for FD mode solver result

%clc; clear variables; close all;
close all;

%load('FD Solver Result.mat');

figure;
plot(real(RetVal.beta) - real(RetVal.beta(1)), 'rx', 'MarkerSize', 8);
xlabel('Mode No.');
ylabel('\Delta\beta');

saveas(gcf, 'Real Prop Constant.png');

figure;
plot(imag(RetVal.beta)*8.68, 'rx', 'MarkerSize', 8);
xlabel('Mode No.');
ylabel('Loss (dB/m)');
ylim([-100 0]);

saveas(gcf, 'Imag Prop Constant.png');