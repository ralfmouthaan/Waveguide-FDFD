% Ralf Mouthaan
% University of Cambridge
% December 2021
%
% Script to plot propagation constant for FD mode solver result

clc; clear variables; close all;

load('FD Solver Result.mat');

figure;
plot(real(RetVal.beta) - real(RetVal.beta(1)), 'rx', 'MarkerSize', 8);
xlabel('Mode No.');
ylabel('\Delta\beta');