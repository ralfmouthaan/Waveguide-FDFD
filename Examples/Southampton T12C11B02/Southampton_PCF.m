% Ralf Mouthaan
% University of Cambridge
% June 2020
%
% Script to run FDFD model on new Erlangen fibre.

clc; clear variables; close all;

%% Set up problem

lambda = 450e-9;
k0 = 2*pi/lambda;
beta = k0*1.3575; % Propagation constant will be close to that of free space.
Nx = 1300;
NoModes = 30;

[x, n] = GenerateFibreProfile(Nx);
dx = x(2) - x(1);

%% Show refractive index profile

figure;
imagesc(x*1e6, x*1e6, n);
axis square;
xlabel('\mum');
ylabel('\mum');
hold on

%% Call FD solver

tic
RetVal = ModeSolverFD(dx, n, lambda, beta, NoModes);
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

    n_silica = 1.45;
    n_pentane = 1.3575;
    r2 = 13.35e-6;
    r3 = 10.55e-6/2;
    s = 350e-9;
    r1 = 47.3e-6/2;
    a = 10.55e-6/2; % Major axis of ellipse
    b = 10.55e-6/2; % Minor axis of ellipse
    
    x = linspace(-r1*1.2,r1*1.2,Nx);
    [x_mesh, y_mesh] = meshgrid(x, x.');
    r_mesh = sqrt(x_mesh.^2 + y_mesh.^2);

    n = ones(Nx, Nx);
    n = n*n_silica;
    n(r_mesh < r1) = n_pentane;
    
    NoCapillaries = 7;
    for ii = 1:NoCapillaries
        alpha = 2*pi/NoCapillaries*ii;
        mux = cos(alpha)*(s+r2+r3);
        muy = sin(alpha)*(s+r2+r3);
        n(...
            (cos(alpha)*(x_mesh - mux) + sin(alpha)*(y_mesh - muy)).^2/(a)^2 + ...
            (sin(alpha)*(x_mesh - mux) - cos(alpha)*(y_mesh - muy)).^2/b^2 <= 1) = n_silica;
        n(...
            (cos(alpha)*(x_mesh - mux) + sin(alpha)*(y_mesh - muy)).^2/(a-s)^2 + ...
            (sin(alpha)*(x_mesh - mux) - cos(alpha)*(y_mesh - muy)).^2/(b-s)^2 <= 1) = n_pentane;
    end
    
    n(r_mesh > r1) = n_silica;

end
