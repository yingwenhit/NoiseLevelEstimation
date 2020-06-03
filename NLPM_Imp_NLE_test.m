clear all
close all

% Input images
name_img = 'mdn4_w';%''
std_n = 10;
load (['NLE_image/' name_img '-In-' num2str(std_n) ])

% Parameters
iter   = 1;%25 185;%20 110;%15 80;%10
K      = 20; 
sigmaG = 1.5; 
tau    = 19; % Parameters in NLPM


% [I,PSNR] = NLPM(I, In, iter, K, sigmaG, tau);
alpha  = 0.8; % alpha<1
beta   = 1e-3; % Parameters in NLE algorithm (see paper)


time = cputime;
[NLE,  u_pde, O_flat, Count] = NLPM_Imp_NLE(In, iter, K, sigmaG, tau, alpha, beta);
time = cputime - time;
dif_noise = std_n - NLE;
disp(['Difference = ' num2str(std_n - NLE) ' Time:' num2str(time) 's'])


% save(['NLE_image/'  name_img '-IMPGS-NLE-' num2str(std_n)])