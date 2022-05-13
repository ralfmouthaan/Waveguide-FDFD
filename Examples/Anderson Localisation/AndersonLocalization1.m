% Ralf Mouthaan
% University of Cambridge
% June 2020
%
% Foray into looking at Anderson localisation using FDFD code.

clc; clear variables; close all;

%% Set up problem

lambda = 633e-9;
k0 = 2*pi/lambda;
beta = k0*1.5; % Propagation constant will be close to that of free space.
Nx = 500;
NoModes = 10;

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

%save('FD Solver Result.mat', 'RetVal', '-v7.3');

function [x,n] = GenerateFibreProfile(Nx)

    x = linspace(-15e-6,15e-6, Nx);
    [x_mesh, y_mesh] = meshgrid(x, x.');
    
    n = ones(Nx);
    idx = rand(Nx);
    n(idx > 0.5) = 1.5;
    
    n(x_mesh>14e-6)= 1;
    n(x_mesh<-14e-6)= 1;
    n(y_mesh>14e-6)= 1;
    n(y_mesh<-14e-6)= 1;
    

end
