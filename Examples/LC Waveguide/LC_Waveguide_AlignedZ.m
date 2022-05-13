% Ralf Mouthaan
% University of Cambridge
% June 2020
%
% Script to run FDFD model on new Erlangen fibre.

clc; clear variables; close all;

%% Set up problem

lambda = 1550e-9;
beta0 = 2*pi*1.5/lambda; % Propagation constant will be close to that of free space.
Nx = 1000;
NoModes = 6;

[x, nx, ny, nz] = GenerateFibreProfile(Nx);
dx = x(2) - x(1);

%% Show refractive index profile

figure;
imagesc(x*1e6, x*1e6, nz);
axis square;
xlabel('\mum');
ylabel('\mum');
hold on

%% Call FD solver

tic
RetVal = ModeSolverFD_Anisotropic(dx, nx, ny, nz, lambda, beta0, NoModes);
toc

%% Plot modes

figure;
for i = 1:6
    
    Ex = RetVal.Ex{i};
    Ey = RetVal.Ey{i};
    Nx = RetVal.Nx;
    x1 = linspace(-1, 1, Nx);
    x2 = linspace(-1, 1, 50);
    Ex = interp2(x1, x1.', Ex, x2, x2.');
    Ey = interp2(x1, x1.', Ey, x2, x2.');

    subplot(3,2,i);
    quiver(x2, x2.', real(Ex), real(Ey));
    hold on
    contour(x1, x1.', RetVal.Eabs{i}/max(max(RetVal.Eabs{i})));
    xlim([-1 1]); ylim([-1 1]);
    xticks(''); yticks('');
    title([num2str(i) '; \beta  = ' num2str(RetVal.beta(i))]);
    
end

fprintf('Mode 1 Effective n = %0.6f\n', RetVal.beta(1)*lambda/2/pi);
fprintf('Mode 2 Effective n = %0.6f\n', RetVal.beta(2)*lambda/2/pi);

function [x, nx, ny, nz] = GenerateFibreProfile(Nx)

    ng = 1.45;
    n0 = 1.5292;
    ne = 1.7072;
    w = 3e-6;
    
    x = linspace(-1.5*w,1.5*w,Nx);
    [x_mesh, y_mesh] = meshgrid(x, x.');
    
    n = ones(Nx, Nx);
    n(y_mesh < w/2) = ng;
    nx = n; ny = n; nz = n;
    nx(abs(x_mesh) <= w/2 & abs(y_mesh) <= w/2) = n0;
    ny(abs(x_mesh) <= w/2 & abs(y_mesh) <= w/2) = n0;
    nz(abs(x_mesh) <= w/2 & abs(y_mesh) <= w/2) = ne;

end