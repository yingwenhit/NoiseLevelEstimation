clear all
close all

addpath('Other Methods NLE/PCANLE_20120521/matlab/')
addpath('Other Methods NLE/NoiseLevel20150203/');

Path = 'NLE_image/';        % ?????
Dir  = dir([Path '*IMPGS-NLE*']); % ????jpg????
nlpm_dif = [];
pca_dif = [];
nlpm_var = [];
pca_var = [];
for i = 1:length(Dir)          % ???????????????
    
%     ImgName = Dir(i).name;
    load([Path Dir(i).name]); %??????
    
    Bch = 100;
    step = 50;
    
    [NLPM_NLE_mini, SIB_NLE_mini, PCA_NLE_mini, NLE_Truth_mini, u_pde, O_flat, Count] ... 
        = NLPM_miniBatch_Imp_NLE_v3(I, In, iter, K, sigmaG, tau, alpha, beta, Bch, step);
    
    nlpm = mean(NLPM_NLE_mini);
    sib = mean(SIB_NLE_mini);
    pca = mean(PCA_NLE_mini);
    truth = mean(NLE_Truth_mini);
    
    nlpm_var = var(NLPM_NLE_mini);
    sib_var = var(SIB_NLE_mini);
    pca_var = var(PCA_NLE_mini);
    truth_var = var(NLE_Truth_mini);
    
    disp([Dir(i).name...
        ' | nlpm: ' num2str(nlpm) ' ' num2str(nlpm_var)...
        ' | pca: ' num2str(pca) ' ' num2str(pca_var)...
        ' | sib: ' num2str(sib) ' ' num2str(sib_var)...
        ' | truth: ' num2str(truth) ' ' num2str(truth_var)]);
    [~, sz] = size(NLE_mini);
%     figure,plot(1:sz, NLPM_NLE_mini, 1:sz, SIB_NLE_mini, 1:sz,...
%         PCA_NLE_mini, 1:sz, NLE_Truth_mini)
%     legend('nlpm','sib','pca','truth')
%     
%     sib_mini = mean(SIB_NLE_mini);
%     sib_mini_dif = std_n-sib_mini;
%     disp([Dir(i).name ' pca: ' num2str(pca_nle) ' pca_dif: ' num2str(pca_dif)...
%          ' sib: ' num2str(sib_nle) ' sib_dif: ' num2str(sib_dif)...
%          ' sib_mini: ' num2str(sib_mini) ' sib_mini_dif: ' num2str(sib_mini_dif)]);
     
     %     pca_nle = sqrt( PCANoiseLevelEstimator(In) );
%     pca_dif = std_n - pca_nle;
%     
%     
%     tic;
%     [nlevel ~] = NoiseLevel(In);
%     t=toc;
%     sib_nle = nlevel;
%     sib_dif = std_n - sib_nle;

%     time = cputime;
%     [NLE_mini, PCA_NLE_mini, u_pde, O_flat, Count] ...
%         = NLPM_miniBatch_Imp_NLE_v2(I, In, iter, K, sigmaG, tau, alpha, beta, Bch, step);
%     time = cputime - time;
%     nle = mean(NLE_mini);
%     nle_var = var(NLE_mini);
%     pca = mean(PCA_NLE_mini);
%     pcavar = var(PCA_NLE_mini);
%     
%     nlpm_dif = [nlpm_dif, std_n-nle];
%     nlpm_var = [nlpm_var, nle_var];
%     
%     pca_dif = [pca_dif, std_n-pca];
%     pca_var = [pca_var, pcavar];
%     
%     disp([Dir(i).name ' nlpm: ' num2str(NLE) ' nlpm_dif: ' num2str(dif_noise) ...
%         ' nlpm_mini: ' num2str(nle) ' nlpm_mini_dif: ' num2str(std_n-nle) ...
%         ' pca_mini: ' num2str(pca) ' pca_mini_dif: ' num2str(std_n-pca)])
% %     dif_noise_ = [dif_noise_; std_n - NLE];
% %     disp(['Difference = ' num2str(std_n - NLE) ' Time:' num2str(time) 's'])
%     
    
%     save(['NLE_image/'  name_img '-IMPGS-miniBatch-NLE-' num2str(std_n)])
    

end
% mean(abs(dif_noise_))

