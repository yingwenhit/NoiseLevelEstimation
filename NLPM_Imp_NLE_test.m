clear all
close all

% Input images
name_img = 'cameraman';%''
std_n = 20;
load (['NLE_image/' name_img '-In-' num2str(std_n) ])

% Parameters
iter   = 1;%25 185;%20 110;%15 80;%10
K      = 5; 
sigmaG = 2; 
tau    = 30; % Parameters in NLPM


% [I,PSNR] = NLPM(I, In, iter, K, sigmaG, tau);
alpha  = 0.6; % alpha<1
beta   = 1e-3; % Parameters in NLE algorithm (see paper)


time = cputime;
[NLE,  u_pde, O_flat, Count] = NLPM_Imp_NLE(In, iter, K, sigmaG, tau, alpha, beta);
time = cputime - time;
dif_noise = std_n - NLE;
disp(['Difference = ' num2str(std_n - NLE) ' Time:' num2str(time) 's'])

Bch = 100;
step = 50;
[NLE_mini, PCA_NLE, u_pde, O_flat, Count] = NLPM_miniBatch_Imp_NLE_v2(I, In, iter, K, sigmaG, tau, alpha, beta, Bch, step);
nle_mini = mean(NLE_mini);
dif_noise_mini = std_n - nle_mini;
disp(['miniDifference = ' num2str(dif_noise_mini)])

% save(['NLE_image/'  name_img '-IMPGS-NLE-' num2str(std_n)])