% Ralf Mouthaan
% University of Cambridge
% June 2020
%
% Script to run FDFD model on new Erlangen fibre.

clc; clear variables; close all;

%% Set up problem geometry

Nx = 1250;
NoModes = 12;

[x, n] = GenerateFibreProfile(Nx);
dx = x(2) - x(1);

%% Show refractive index profile

figure;
imagesc(x*1e6, x*1e6, n);
axis square;
xlabel('\mum');
ylabel('\mum');
hold on
pause(1);

%% Call FD solver

fid = fopen('Wavelength Sweep - Prop Constants.txt','w+');

for lambda = 400e-9:20e-9:1000e-9
    k0 = 2*pi/lambda;
    fprintf('Wavelength = %0.1fnm\n', lambda*1e9);
    RetVal = ModeSolverFD(dx, n, lambda, k0, NoModes);
    fprintf(fid, 'Wavelength = %0.1fnm\n', lambda*1e9);
    fprintf(fid, 'Mode 1-2: %0.5f + %0.5fi\n', real(mean(RetVal.beta(1:2))), imag(mean(RetVal.beta(1:2))));
    fprintf(fid, 'Mode 3-6: %0.5f + %0.5fi\n', real(mean(RetVal.beta(3:6))), imag(mean(RetVal.beta(3:6))));
    fprintf(fid, 'Mode 7-10: %0.5f + %0.5fi\n', real(mean(RetVal.beta(7:10))), imag(mean(RetVal.beta(7:10))));
    fprintf(fid, 'Mode 11-12: %0.5f + %0.5fi\n\n', real(mean(RetVal.beta(11:12))), imag(mean(RetVal.beta(11:12))));
end

fclose(fid);

function [x,n] = GenerateFibreProfile(Nx)

    n_silica = 1.45;
    r2 = 29.4e-6/2;
    r3 = 9.96e-6/2;
    s = 200e-9;
    r1 = r2 + 2*r3 + 2*s;
    a = 10.5e-6/2; % Major axis of ellipse
    b = 9.73e-6/2; % Minor axis of ellipse
    
    x = linspace(-r1*1.2,r1*1.2,Nx);
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
