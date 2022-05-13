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

%% Parameters for interpolation

Nx1 = size(RetVal.Ex{1}, 1);
Nx2 = 30;
x1 = ((-Nx1/2+0.5):(Nx1/2-0.5))*RetVal.dx;
x2 = linspace(min(x1), max(x1), Nx2);
[xmesh1, ymesh1] = meshgrid(x1, x1.');
[xmesh2,ymesh2] = meshgrid(x2,x2.');

arrRotAngle = linspace(0, 2*pi, 16);

figure;
for i = 1:length(arrRotAngle)
    rotangle = arrRotAngle(i);
    Ex = RetVal.Ex{13} + exp(1i*rotangle)*RetVal.Ex{14};
    Ey = RetVal.Ey{13} + exp(1i*rotangle)*RetVal.Ey{14};
    Eabs = sqrt(abs(Ex).^2 + abs(Ey).^2);

    subplot(4,4,i);

    Ex = interp2(xmesh1, ymesh1, Ex, xmesh2, ymesh2);
    Ey = interp2(xmesh1, ymesh1, Ey, xmesh2, ymesh2);

    contour(x1*1e6, x1.'*1e6, Eabs/max(max(Eabs)));
    axis square;
    hold on
    quiver(x2*1e6, x2.'*1e6, real(Ex), real(Ey));
    xlim([-25 25]);
    ylim([-25 25]);
    xticks('');
    yticks('');
    title(['RotAngle = ' num2str(rotangle)]);
end

%%

rotangle = 4.6077;
Ex = RetVal.Ex{3} + exp(1i*rotangle)*RetVal.Ex{5};
Ey = RetVal.Ey{3} + exp(1i*rotangle)*RetVal.Ey{5};
Eabs = sqrt(abs(Ex).^2 + abs(Ey).^2);

figure;

Ex = interp2(xmesh1, ymesh1, Ex, xmesh2, ymesh2);
Ey = interp2(xmesh1, ymesh1, Ey, xmesh2, ymesh2);

contour(x1*1e6, x1.'*1e6, Eabs/max(max(Eabs)));
axis square;
hold on
quiver(x2*1e6, x2.'*1e6, real(Ex), real(Ey));
xlim([-25 25]);
ylim([-25 25]);
xticks('');
yticks('');
title(['RotAngle = ' num2str(rotangle)]);
    