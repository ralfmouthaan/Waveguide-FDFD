% Ralf Mouthaan
% University of Cambridge
% June 2020
%
% Script to run FDFD model on new Erlangen fibre.

clc; clear variables; close all;

%% Set up problem

lambda = 1500e-9;
beta0 = 2*pi*1.55/lambda; % Propagation constant will be close to that of free space.
Nx = 1000;
NoModes = 9;

[x, nx, ny, nz] = GenerateFibreProfile(Nx);
dx = x(2) - x(1);

%% Show refractive index profile

figure('position', [400 400 400 400]);
imagesc(x*1e6, x*1e6, nz);
axis square;
xlabel('\mum');
ylabel('\mum');

%% Call FD solver

tic
RetVal = ModeSolverFD_Anisotropic(dx, nx, ny, nz, lambda, beta0, NoModes);
toc

%% Plot modes

for i = 1:6
    
    Ex = RetVal.Ex{i};
    Ey = RetVal.Ey{i};
    Nx = RetVal.Nx;
    x1 = (-Nx/2+1/2:Nx/2-1/2)*dx;
    x2 = linspace(min(x1), max(x1), 50);
    Ex = interp2(x1, x1.', Ex, x2, x2.');
    Ey = interp2(x1, x1.', Ey, x2, x2.');

    figure
    quiver(x2*1e6, x2.'*1e6, real(Ex), real(Ey));
    hold on
    contour(x1*1e6, x1.'*1e6, RetVal.Eabs{i}/max(max(RetVal.Eabs{i})));
    xlim([-4 4]); ylim([-4 4]);
    %title(['\beta  = ' num2str(RetVal.beta(i))]);
    axis square;
    xlabel('\mum'); ylabel('\mum');
    xticks(-4:2:4); yticks(-4:2:4);
    
end

fprintf('Mode 1 Effective n = %0.6f\n', RetVal.beta(1)*lambda/2/pi);
fprintf('Mode 2 Effective n = %0.6f\n', RetVal.beta(2)*lambda/2/pi);

function [x, nx, ny, nz] = GenerateFibreProfile(Nx)

    n0 = 1.45;
    ne = 1.55;
    w = 6e-6;
    theta=0; % angle of LC with x axis
    
    x = linspace(-1.5*w,1.5*w,Nx);
    [x_mesh, y_mesh] = meshgrid(x, x.');
    r_mesh = sqrt(x_mesh.^2 + y_mesh.^2);
    
    n = ones(Nx, Nx);
    nx = n; ny = n; nz = n;
    nx(abs(r_mesh) <= w/2) = sqrt(n0^2*cos(theta)^2 + ne^2*sin(theta)^2);
    ny(abs(r_mesh) <= w/2) = sqrt(n0^2*sin(theta)^2 + ne^2*cos(theta)^2);
    nz(abs(r_mesh) <= w/2) = n0;

end