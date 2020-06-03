% WRITE THE NLE IMP GS RESULTS
clear; close all;

imgDir  = dir('*IMPGS-NLE*.mat'); 
name_f = 'IMPGS_NLE_results_parameter.txt';
fp = fopen(name_f,'w');

fprintf(fp, 'NLE  dif_noise name  iter  K   sigmaG  tau   alpha   beta \n');
% disp([ 'NLE  dif_noise name  iter  K   sigmaG  tau   alpha   beta']);

for i = 1:length(imgDir)
    
    load([imgDir(i).name]); 
    fprintf(fp, '%s, %s, %s, %s, %s, %s, %s, %s, %s\n', num2str(NLE),num2str(dif_noise),[imgDir(i).name],num2str(iter),num2str(K),num2str(sigmaG),...
         num2str(tau),num2str(alpha),num2str(beta));
%     disp([num2str(NLE) ' ' num2str(dif_noise) ' '  [imgDir(i).name] ' ' num2str(iter) ' ' num2str(K) ' ' num2str(sigmaG) ' '...
%          num2str(tau) ' ' num2str(alpha) ' ' num2str(beta)]);

end
fclose(fp);