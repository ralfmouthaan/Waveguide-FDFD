% Ralf Mouthaan
% University of Cambridge
% June 2020
%
% Script to run FDFD model on new Erlangen fibre.

clc; clear variables; %close all;

%% Set up problem

lambda = 633e-9; % Wavelength in nm
k0 = 2*pi/lambda; % Wavenumber
beta = k0; % Fundamental mode propagation constant will be close to that of free space.
Nx = 1000; % Size of problem space. Assumed to be square, i.e. Nx x Nx
NoModes = 30; % Number of modes we wish to find.

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
    r2 = 29.4e-6/2;
    r3 = 9.96e-6/2;
    s = 200e-9;
    r1 = r2 + 2*r3 + 2*s;
    a = 10.5e-6/2; % Major axis of ellipse
    b = 9.73e-6/2; % Minor axis of ellipse
    
    x = linspace(-r1*1.5,r1*1.5,Nx);
    [x_mesh, y_mesh] = meshgrid(x, x.');
    r_mesh = sqrt(x_mesh.^2 + y_mesh.^2);

    n = ones(Nx, Nx);
    n = n*n_silica;
    n(r_mesh < r1) = 1;
    
    NoCapillaries = 8;
    for ii = 1:NoCapillaries
        alpha = 2*pi/NoCapillaries*ii;
        mux = cos(alpha)*(s+r2+r3*1.03);
        muy = sin(alpha)*(s+r2+r3*1.03);
        n(...
            (cos(alpha)*(x_mesh - mux) + sin(alpha)*(y_mesh - muy)).^2/(a*1.03)^2 + ...
            (sin(alpha)*(x_mesh - mux) - cos(alpha)*(y_mesh - muy)).^2/b^2 <= 1) = n_silica;
        n(...
            (cos(alpha)*(x_mesh - mux) + sin(alpha)*(y_mesh - muy)).^2/(a*1.03-s)^2 + ...
            (sin(alpha)*(x_mesh - mux) - cos(alpha)*(y_mesh - muy)).^2/(b-s)^2 <= 1) = 1;
    end
    
    n(r_mesh > r1) = n_silica;

end
