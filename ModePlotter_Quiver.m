% Ralf Mouthaan
% University of Cambridge
% August 2021
%
% Script to superimpose Southampton Fibre modes to find linearly polarised
% ones.

clc; close all;
clearvars -except RetVal;

if exist('RetVal', 'var') == 0
    load('FD Solver Result.mat');
end

%% Parameters

Nx1 = size(RetVal.Ex{1}, 1);
Nx2 = 30;

%% Propagation Constants

figure;
subplot(1,2,1);
plot(real(RetVal.beta)/RetVal.k0, 'rx');
xlim([1 20]);
xlabel('Mode Number');
ylabel('Effective n');
subplot(1,2,2);
plot(imag(RetVal.beta)*20/log(10), 'rx');
xlim([1 20]);
xlabel('Mode Number');
ylabel('Loss (dB/m)');

%% Plot Raw Modes

% Interpolation
x1 = ((-Nx1/2+0.5):(Nx1/2-0.5))*RetVal.dx;
x2 = linspace(min(x1), max(x1), Nx2);
[xmesh1, ymesh1] = meshgrid(x1, x1.');
[xmesh2,ymesh2] = meshgrid(x2,x2.');

figure;
for i = 1:16
    
    subplot(4,4,i);

    Ex = interp2(xmesh1, ymesh1, RetVal.Ex{i}, xmesh2, ymesh2);
    Ey = interp2(xmesh1, ymesh1, RetVal.Ey{i}, xmesh2, ymesh2);
    Eabs = RetVal.Eabs{i};

    contour(x1*1e6, x1.'*1e6, Eabs/max(max(Eabs)));
    axis square;
    hold on
    quiver(x2*1e6, x2.'*1e6, real(Ex), real(Ey));
    xlim([-25 25]);
    ylim([-25 25]);
    xticks('');
    yticks('');
    
end