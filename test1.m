%% Experiment 1
% We generate instances of the blind deconvolution problem and initialize at
% random points in a hypercube. Then we run polyak's method see how many instances recover a point in the
% solution set.
clc; clear all; close all;
d1 = 100;
d2 = 50;
Cs = 1:8;
nus = 2.^(2:10);
num_rep = 10;
res_error = zeros(length(Cs), length(nus),2);
res_count = zeros(length(Cs), length(nus),2);
res_iter = zeros(length(Cs), length(nus), num_rep,2);


for ii = 1:length(Cs)
    m = Cs(ii)*(d1+d2);
    for jj = 1:length(nus)
        nu = nus(jj);
        fprintf('Constant %d, nu %7.2e \n', Cs(ii), nu);
        for kk = 1:num_rep
            fprintf('\t Instance %d \n', kk);
            wb = eye(d1,1);
            xb = eye(d2,1);
            
            L = randn(m,d1);
            R = randn(m,d2);
            y = (L*wb).*(R*xb);
            
            % Random init on hypercube
            w0 = 2*nu*(rand(d1,1)-0.5);
            x0 = 2*nu*(rand(d2,1)-0.5);
            [w,x, f] = BD_polyak(L,R,y,w0,x0);
            res_error(ii,jj,1) = res_error(ii,jj,1)+ (norm(w*x'-wb*xb','fro')/norm(wb*xb','fro'))/num_rep;
            res_count(ii,jj,1) = res_count(ii,jj,1)+ (norm(w*x'-wb*xb','fro')/norm(wb*xb','fro') <= 1e-6)/num_rep;
            res_iter(ii,jj,kk,1) = length(f);
            
            % Gaussian init
            w0 = nu*(randn(d1,1)/sqrt(d1));
            x0 = nu*(randn(d2,1)/sqrt(d2));
            [w,x, f] = BD_polyak(L,R,y,w0,x0);
            res_error(ii,jj,2) = res_error(ii,jj,2)+ (norm(w*x'-wb*xb','fro')/norm(wb*xb','fro'))/num_rep;
            res_count(ii,jj,2) = res_count(ii,jj,2)+ (norm(w*x'-wb*xb','fro')/norm(wb*xb','fro') <= 1e-6)/num_rep;
            res_iter(ii,jj,kk,2) = length(f);
            save results.mat;
        end
    end
end
 


% clc; clear all; close all;
% fprintf('RUNNING FOR 200 \n')
% d1 = 200;
% d2 = 100;
% Cs = 1:8;
% nus = 2.^(2:10);
% num_rep = 10;
% res_error = zeros(length(Cs), length(nus),2);
% res_count = zeros(length(Cs), length(nus),2);
% res_iter = zeros(length(Cs), length(nus), num_rep,2);
% 
% 
% for ii = 1:length(Cs)
%     m = Cs(ii)*(d1+d2);
%     for jj = 1:length(nus)
%         nu = nus(jj);
%         fprintf('Constant %d, nu %7.2e \n', Cs(ii), nu);
%         for kk = 1:num_rep
%             fprintf('\t Instance %d \n', kk);
%             wb = eye(d1,1);
%             xb = eye(d2,1);
%             
%             L = randn(m,d1);
%             R = randn(m,d2);
%             y = (L*wb).*(R*xb);
%             
%             % Random init on hypercube
%             w0 = 2*nu*(rand(d1,1)-0.5);
%             x0 = 2*nu*(rand(d2,1)-0.5);
%             [w,x, f] = BD_polyak(L,R,y,w0,x0);
%             res_error(ii,jj,1) = res_error(ii,jj,1)+ (norm(w*x'-wb*xb','fro')/norm(wb*xb','fro'))/num_rep;
%             res_count(ii,jj,1) = res_count(ii,jj,1)+ (norm(w*x'-wb*xb','fro')/norm(wb*xb','fro') <= 1e-6)/num_rep;
%             res_iter(ii,jj,kk,1) = length(f);
%             
%             % Gaussian init
%             w0 = nu*(randn(d1,1)/sqrt(d1));
%             x0 = nu*(randn(d2,1)/sqrt(d2));
%             [w,x, f] = BD_polyak(L,R,y,w0,x0);
%             res_error(ii,jj,2) = res_error(ii,jj,2)+ (norm(w*x'-wb*xb','fro')/norm(wb*xb','fro'))/num_rep;
%             res_count(ii,jj,2) = res_count(ii,jj,2)+ (norm(w*x'-wb*xb','fro')/norm(wb*xb','fro') <= 1e-6)/num_rep;
%             res_iter(ii,jj,kk,2) = length(f);
%             save results200pow.mat;
%         end
%     end
% end
