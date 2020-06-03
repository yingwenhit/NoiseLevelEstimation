clear;close all

imgDir  = dir('../../NLE_image/*-In-*.mat');

name_f = 'NoiseLevel2015_results.txt';
fp = fopen(name_f,'w');


for i=1:length(imgDir)
    load(['../../NLE_image/' imgDir(i).name]);
    noise = In;
    %  noise = img + randn(size(img)) * level(i);
    tic;
    [nlevel ~] = NoiseLevel(noise);
    t=toc;
    
    fprintf(fp, '%f, %s \n', nlevel, [imgDir(i).name]);
    
    
end
fclose(fp);