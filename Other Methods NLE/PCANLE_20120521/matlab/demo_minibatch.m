%%%%%%%%%%%%%%%%% experiments test, 2020/1/7 %%%%%%%%%%%%%%%%%%%%%%%
clear all;close all
% imgDir  = dir('../../../NLE_image/*-In-*.mat');
imgDir = dir('../../../Lena-10-minitest/*.mat');

% name_f = 'PCANLE2012_results.txt';
% fp = fopen(name_f,'w');
sigma = [];

for i = 1:length(imgDir)
    
    load(['../../../Lena-10-minitest/' imgDir(i).name]);
    x = f;
    
    if size(x,3) ~= 1
        error( 'Only grayscale 2D images can be processed.' );
    end
    
    sigma_ = sqrt( PCANoiseLevelEstimator(x) );
    sigma = [sigma, sigma_];
%     disp('sigma');
%     fprintf(fp, '%f, %s \n', sigma, [imgDir(i).name]);
end

% fclose(fp);