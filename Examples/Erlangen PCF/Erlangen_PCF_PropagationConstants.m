% Ralf Mouthaan
% University of Cambridge
%  June 2021
%
% Plotting out the propagation constants.

clc; close all;

if exist('RetVal', 'var') == 0
    load('FD Solver Result - 3rd Order PML.mat');
end

%% Real component

figure;
plot(1:2,real(RetVal.beta(1:2))/RetVal.k0, 'x', 'MarkerSize', 8);
hold on
plot(3:6,real(RetVal.beta(3:6))/RetVal.k0, 'o', 'MarkerSize', 8);
plot(7:10,real(RetVal.beta(7:10))/RetVal.k0, '*', 'MarkerSize', 8);
plot(11:12,real(RetVal.beta(11:12))/RetVal.k0, 's', 'MarkerSize', 8);
plot(13:16,real(RetVal.beta(13:16))/RetVal.k0, '+', 'MarkerSize', 8);
plot(17:20,real(RetVal.beta(17:20))/RetVal.k0, 'd', 'MarkerSize', 8);
plot(21:22,real(RetVal.beta(21:22))/RetVal.k0, '^', 'MarkerSize', 8);
plot(23:30,real(RetVal.beta(23:30))/RetVal.k0, 'p', 'MarkerSize', 8);
xlabel('Mode No.');
ylabel('Effective n');
ylim([0.9988 1]);
set(gca, 'FontSize', 12);

legend('LP_{0,1}', 'LP_{1,1}', 'LP_{2,1}', 'LP_{0,2}', 'LP_{3,1}', 'LP_{1,2}', 'LP_{4,1}', 'Non-Guided');

%% Imaginary component

figure
semilogy(1:2,-imag(RetVal.beta(1:2))*20/log(10), 'x', 'MarkerSize', 8);
hold on
semilogy(3:6,-imag(RetVal.beta(3:6))*20/log(10), 'o', 'MarkerSize', 8);
semilogy(7:10,-imag(RetVal.beta(7:10))*20/log(10),  '*', 'MarkerSize', 8);
semilogy(11:12,-imag(RetVal.beta(11:12))*20/log(10),  's', 'MarkerSize', 8);
semilogy(13:16,-imag(RetVal.beta(13:16))*20/log(10),  '+', 'MarkerSize', 8);
semilogy(17:20,-imag(RetVal.beta(17:20))*20/log(10), 'd', 'MarkerSize', 8);
semilogy(21:22,-imag(RetVal.beta(21:22))*20/log(10),  '^', 'MarkerSize', 8);
semilogy(23:30,-imag(RetVal.beta(23:30))*20/log(10),  'p', 'MarkerSize', 8);
xlabel('Mode No.');
ylabel('Modal Loss (dB/m)');
set(gca, 'FontSize', 12);
yticks([1 10 100 800]);

legend('LP_{0,1}', 'LP_{1,1}', 'LP_{2,1}', 'LP_{0,2}', 'LP_{3,1}', 'LP_{1,2}', 'LP_{4,1}', 'Non-Guided');