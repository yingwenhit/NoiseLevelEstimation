% WRITE THE NLE RESULTS
clear; close all;

imgDir  = dir('*NLE*.mat'); 
name_f = 'NLE_results_parameter.txt';
fp = fopen(name_f,'w');

for i = 1:length(imgDir)
    
    load([imgDir(i).name]); 
%     fprintf(fp, '%f, %f, %s\n', NLE, dif_noise, [imgDir(i).name]);
%     disp([num2str(NLE) ' ' num2str(dif_noise) ' '  [imgDir(i).name] ' ' num2str(iter) ' ' num2str(K) ' ' num2str(sigmaG) ' '...
%          num2str(tau) ' ' num2str(alpha) ' ' num2str(beta)]);
    disp([num2str(beta)])
%     iter   = 560;%25 185;%20 110;%15 80;%10
% K      = 2; 
% sigmaG = 1; 
% tau    = 0.2; % Parameters in NLPM
% 
% 
% % [I,PSNR] = NLPM(I, In, iter, K, sigmaG, tau);
% alpha  = 0.2; 
% beta   = 1e-3;
end
fclose(fp);