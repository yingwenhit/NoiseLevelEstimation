clear all
close all

% Input images
name_img = 'lena';%''
std_n = 10;
load (['NLE_image/' name_img '-In-' num2str(std_n) ])

% Parameters
iter   = 1;%25 185;%20 110;%15 80;%10
K      = 8;
sigmaG = 1.5;
tau    = 13; % Parameters in NLPM

% [I,PSNR] = NLPM(I, In, iter, K, sigmaG, tau);
alpha  = 0.8; % alpha<1
beta   = 1e-3; % Parameters in NLE algorithm (see paper)


Bch = 60;
step = 50;
time = cputime;
[NLE,  u_pde, O_flat, Count] = NLPM_miniBatch_Imp_NLE_v2(I, In, iter, K, sigmaG, tau, alpha, beta, Bch, step);
time = cputime - time;
disp(['Difference = ' num2str(std_n - NLE) ' Time:' num2str(time) 's'])


%     save(['NLE_image/'  name_img '-IMPGS-miniBatch-NLE-' num2str(std_n)])


