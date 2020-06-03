%% Test for same parameters
clear all
close all

%% Parameters
iter   = 400;%25 185;%20 110;%15 80;%10
K      = 1;
sigmaG = 1;
tau    = 0.2; % Parameters in NLPM
alpha  = 0.2;
beta   = 1e-3; % Parameters in NLE algorithm (see paper)

% Input images
Name = {'mdn4_w', 'cameraman', 'fingerprint1', 'lena'};
Std_n = {10, 15, 20, 25};
for i=1:4
    for j=1:4
        name_img = Name{i};
        std_n = Std_n{j};
        load (['NLE_image/' name_img '-In-' num2str(std_n) ])
        
        [NLE,  u_pde, O_flat, Count] = NLPM_NLE_C2(In, iter, K, sigmaG, tau, alpha, beta);
        dif_noise = std_n - NLE;
%         disp(['Difference = ' num2str(std_n - NLE)])
        disp([name_img ' ' num2str(std_n) ' ' num2str(dif_noise)]);
    end
end




% save(['NLE_image/'  name_img '-NLE-' num2str(std_n)])