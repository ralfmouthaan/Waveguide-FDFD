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

%% Interpolation parameters

Nx1 = size(RetVal.Ex{1}, 1);
Nx2 = 30;
x1 = ((-Nx1/2+0.5):(Nx1/2-0.5))*RetVal.dx;
x2 = linspace(min(x1), max(x1), Nx2);
[xmesh1, ymesh1] = meshgrid(x1, x1.');
[xmesh2,ymesh2] = meshgrid(x2,x2.');

%% Calculate New Modes
Ex_New{1} = RetVal.Ex{1};
Ey_New{1} = RetVal.Ey{1};
Eabs_New{1} = RetVal.Eabs{1};
beta_New(1) = RetVal.beta(1);

Ex_New{2} = RetVal.Ex{5} + exp(1i*2.5133)*RetVal.Ex{6};
Ey_New{2} = RetVal.Ey{5} + exp(1i*2.5133)*RetVal.Ey{6};
Eabs_New{2} = sqrt(abs(Ex_New{2}).^2 + abs(Ey_New{2}).^2);
beta_New(2) = (RetVal.beta(5) + RetVal.beta(6))/2;

Ex_New{3} = RetVal.Ex{3} + exp(1i*4.61)*RetVal.Ex{5};
Ey_New{3} = RetVal.Ey{3} + exp(1i*4.61)*RetVal.Ey{5};
Eabs_New{3} = sqrt(abs(Ex_New{3}).^2 + abs(Ey_New{3}).^2);
beta_New(3) = (RetVal.beta(3) + RetVal.beta(5))/2;

Ex_New{4} = RetVal.Ex{7};
Ey_New{4} = RetVal.Ey{7};
Eabs_New{4} = RetVal.Eabs{7};
beta_New(4) = RetVal.beta(7);

Ex_New{5} = RetVal.Ex{8} + exp(1i*0)*RetVal.Ex{9};
Ey_New{5} = RetVal.Ey{8} + exp(1i*0)*RetVal.Ey{9};
Eabs_New{5} = sqrt(abs(Ex_New{5}).^2 + abs(Ey_New{5}).^2);
beta_New(5) = (RetVal.beta(8) + RetVal.beta(9))/2;

Ex_New{6} = RetVal.Ex{11};
Ey_New{6} = RetVal.Ey{11};
Eabs_New{6} = sqrt(abs(Ex_New{6}).^2 + abs(Ey_New{6}).^2);
beta_New(6) = RetVal.beta(11);

Ex_New{7} = RetVal.Ex{13};
Ey_New{7} = RetVal.Ey{13};
Eabs_New{7} = RetVal.Eabs{13};
beta_New(7) = RetVal.beta(13);

Ex_New{8} = RetVal.Ex{14};
Ey_New{8} = RetVal.Ey{14};
Eabs_New{8} = RetVal.Eabs{14};
beta_New(8) = RetVal.beta(14);

RetVal_New = RetVal;
RetVal_New.Ex = Ex_New;
RetVal_New.Ey = Ey_New;
RetVal_New.Eabs = Eabs_New;
RetVal_New.beta = beta_New;

%% Plot New Modes

figure;
for i = 1:length(RetVal_New.Ex)
    
    if isempty(RetVal_New.Ex{i})
        continue
    end
    subplot(2,4,i);

    Ex = interp2(xmesh1, ymesh1, RetVal_New.Ex{i}, xmesh2, ymesh2);
    Ey = interp2(xmesh1, ymesh1, RetVal_New.Ey{i}, xmesh2, ymesh2);
    Eabs = RetVal_New.Eabs{i};

    contour(x1*1e6, x1.'*1e6, Eabs/max(max(Eabs)));
    axis square;
    hold on
    quiver(x2*1e6, x2.'*1e6, real(Ex), real(Ey));
    xlim([-25 25]);
    ylim([-25 25]);
    xticks('');
    yticks('');
    title(['Mode #' num2str(i)]);
    
end

%% Save new modes
    