
close all;

%load('FD Solver Result.mat');

for i = 1:40
    figure;
    imagesc(x*1e6, x*1e6, RetVal.Eabs{i});
    %title(['\beta = ' num2str(RetVal.beta(i))]);
    axis square;
    xlabel('\mum');
    ylabel('\mum');
    saveas(gcf, ['Mode Profile ' num2str(i) '.png']);
 end