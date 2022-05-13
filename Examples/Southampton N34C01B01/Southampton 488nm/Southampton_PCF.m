% Ralf Mouthaan
% University of Cambridge
% June 2020
%
% Script to run FDFD model on new Erlangen fibre.

clc; clear variables; close all;

%% Set up problem

lambda = 488e-9;
k0 = 2*pi/lambda;
beta = k0*1.3575; % Propagation constant will be close to that of free space.
Nx = 2000;
NoModes = 40;

[x, n] = GenerateFibreProfile(Nx);
dx = x(2) - x(1);

%% Show refractive index profile

figure;
imagesc(x*1e6, x*1e6, n);
axis square;
xlabel('\mum');
ylabel('\mum');

%% Call FD solver

tic
RetVal = ModeSolverFD_LP(dx, n, lambda, beta, NoModes);
toc

%% Plot modes

figure; plot(real(RetVal.beta)/k0);
figure; plot(imag(RetVal.beta)*20/log(10));

for i = 1:NoModes
    figure;
     imagesc(x*1e6, x*1e6, RetVal.Eabs{i});
     title(['\beta = ' num2str(RetVal.beta(i))]);
     axis square;
 end

save('FD Solver Result.mat', 'RetVal', '-v7.3');

function [x,n] = GenerateFibreProfile(Nx)

    % Am assuming this is the first fibre from the Nature paper 
    % https://www.nature.com/articles/s41467-020-19910-7.pdf, as these
    % dimensions best match the dimensions estimated from the SEM +
    % transmission spectrum.

    n_silica = 1.45;
    n_pentane = 1.3575;
    s1 = 750e-9; % large capillary thickness
    s2 = 800e-9; % small capillary thickness
    r = 65e-6/2;
    a = 8.9e-6/2; % Diameter of small capillary
    b = 18.2e-6/2; % Diameter of large capillary
    r1 = r - a + 1e-6; % Offset of small capillary from centre
    r2 = r - b; % Offset of large capillary from centre
    
    x = linspace(-r*1.2,r*1.2,Nx);
    [x_mesh, y_mesh] = meshgrid(x, x.');
    r_mesh = sqrt(x_mesh.^2 + y_mesh.^2);

    n = ones(Nx, Nx);
    n = n*n_silica;
    n(r_mesh < r) = n_pentane;
    
    NoCapillaries = 6;
    
    for ii = 1:NoCapillaries
        alpha = 2*pi/NoCapillaries*ii;
        mux = cos(alpha)*(s1+r2);
        muy = sin(alpha)*(s1+r2);
        n(...
            (cos(alpha)*(x_mesh - mux) + sin(alpha)*(y_mesh - muy)).^2/(b)^2 + ...
            (sin(alpha)*(x_mesh - mux) - cos(alpha)*(y_mesh - muy)).^2/b^2 <= 1) = n_silica;
        n(...
            (cos(alpha)*(x_mesh - mux) + sin(alpha)*(y_mesh - muy)).^2/(b-s1)^2 + ...
            (sin(alpha)*(x_mesh - mux) - cos(alpha)*(y_mesh - muy)).^2/(b-s1)^2 <= 1) = n_pentane;
    end
    
    for ii = 1:NoCapillaries
        alpha = 2*pi/NoCapillaries*ii;
        mux = cos(alpha)*(s2+r1);
        muy = sin(alpha)*(s2+r1);
        n(...
            (cos(alpha)*(x_mesh - mux) + sin(alpha)*(y_mesh - muy)).^2/(a)^2 + ...
            (sin(alpha)*(x_mesh - mux) - cos(alpha)*(y_mesh - muy)).^2/a^2 <= 1) = n_silica;
        n(...
            (cos(alpha)*(x_mesh - mux) + sin(alpha)*(y_mesh - muy)).^2/(a-s2)^2 + ...
            (sin(alpha)*(x_mesh - mux) - cos(alpha)*(y_mesh - muy)).^2/(a-s2)^2 <= 1) = n_pentane;
    end
    
    n(r_mesh > r) = n_silica;

end
