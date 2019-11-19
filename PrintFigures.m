%% Print figures d = 100
clear all; close all; clc;
load('results100pow.mat');
imagesc(res_count(:,:,1))
xticklabels({'2^2','2^3','2^4','2^5','2^6','2^7','2^8','2^9', '2^{10}'})
% xticklabels({'2','4','8','16','32','128','256','512','1024'})
fig = gca;
fig.FontSize = 14; 
xlabel('Scaling \nu','FontSize', 20);
ylabel('$m/(d_1+d_2)$','Interpreter','latex','FontSize', 20);
title('$d_1=100, d_2 =50$', 'Interpreter','latex','FontSize', 20);
colormap gray
export_fig PhaseTransitionCube100.pdf -transparent

imagesc(res_count(:,:,2))
xticklabels({'2^2','2^3','2^4','2^5','2^6','2^7','2^8','2^9', '2^{10}'})
% xticklabels({'2','4','8','16','32','128','256','512','1024'})
fig = gca;
fig.FontSize = 14; 
xlabel('Scaling \nu', 'FontSize', 20);
ylabel('$m/(d_1+d_2)$','Interpreter','latex','FontSize', 20);
title('$d_1=100, d_2 =50$', 'Interpreter','latex','FontSize', 20);
colormap gray
export_fig PhaseTransitionGaussian100.pdf -transparent

%% Print figures d = 200
clear all; close all; clc;
load('results200pow.mat');
imagesc(res_count(:,:,1))
xticklabels({'2^2','2^3','2^4','2^5','2^6','2^7','2^8','2^9', '2^{10}'})
% xticklabels({'2','4','8','16','32','128','256','512','1024'})
fig = gca;
fig.FontSize = 14; 
xlabel('Scaling \nu','FontSize', 20);
ylabel('$m/(d_1+d_2)$','Interpreter','latex','FontSize', 20);
title('$d_1=200, d_2 =100$', 'Interpreter','latex','FontSize', 20);
colormap gray
export_fig PhaseTransitionCube200.pdf -transparent

figure
imagesc(res_count(:,:,2))
xticklabels({'2^2','2^3','2^4','2^5','2^6','2^7','2^8','2^9', '2^{10}'})
% xticklabels({'2','4','8','16','32','128','256','512','1024'})
fig = gca;
fig.FontSize = 14; 
xlabel('Scaling \nu', 'FontSize', 20);
ylabel('$m/(d_1+d_2)$','Interpreter','latex','FontSize', 20);
title('$d_1=200, d_2 =100$', 'Interpreter','latex','FontSize', 20);
colormap gray
export_fig PhaseTransitionGaussian200.pdf -transparent

%% Prints 2D Population
clc; clear all; close all;
grid = linspace(-2,2,500);
z = zeros(length(grid),length(grid));
for ii = 1:length(grid)
    for jj = 1:length(grid)
        w = grid(ii);
        x = grid(jj);
        X = w*x - 1;
        d1 = norm(X);
        d2 = 0;
        [~,E] = ellipke(1-(d2^2/d1^2));
        z(ii,jj) = 2*d1*E/pi;

    end
end
%%
close all;
% colormap winter 
surf(grid, grid, z, 'EdgeColor', 'none');
% set(gca,'xtick', []); 
set(gca, 'xticklabel', []);
% set(gca,'ytick', []); 
set(gca, 'yticklabel', []);
% set(gca,'ztick', []); 
set(gca, 'zticklabel', []);
set(gca, 'box', 'off');
export_fig function2D.pdf -transparent
%%
figure;
 imagesc(z); hold on
% contour(z,20);
set(gca,'xtick', []); 
set(gca, 'xticklabel', []);
set(gca,'ytick', []); 
% axis equal
set(gca, 'yticklabel', []);
export_fig contour2D.pdf -transparent

% colormap winter 
