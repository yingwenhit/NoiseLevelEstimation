clear all
close all

% Input images
name_img = 'lena';
std_n = 25;
load (['NLE_image/' name_img '-In-' num2str(std_n) ])

% Parameters
iter   = 161;%25 185;%20 110;%15 80;%10
K      = 2; 
sigmaG = 1; 
tau    = 0.2; % Parameters in NLPM


% [I,PSNR] = NLPM(I, In, iter, K, sigmaG, tau);
alpha  = 0.2; 
beta   = 1e-3; % Parameters in NLE algorithm (see paper)

time = cputime;
[NLE,  u_pde, O_flat, Count] = NLPM_NLE_C2(In, iter, K, sigmaG, tau, alpha, beta);
time = cputime - time;
dif_noise = std_n - NLE;
disp(['Difference = ' num2str(std_n - NLE) ' Time:' num2str(time) 's'])

% save(['NLE_image/'  name_img '-NLE-' num2str(std_n)])