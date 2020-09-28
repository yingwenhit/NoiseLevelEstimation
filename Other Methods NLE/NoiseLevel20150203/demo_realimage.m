%%%%%%%%%%%%%%%%% experiments test, 2020/1/7 %%%%%%%%%%%%%%%%%%%%%%%
%%% for real image, 2020-06-16
datapath = '../../../NPSID_Results/';
imgDir  = dir(datapath);

name_f = 'SIB_realimage_results.txt';
fp = fopen(name_f,'w');

for i = 1:length(imgDir)
    
    ImgName = imgDir(i).name;
    if (contains(ImgName,'.mat'))
        load([datapath ImgName]);
        
        [sigma1 ~] = NoiseLevel(I1);
        [sigma2 ~] = NoiseLevel(I2);
        [sigma3 ~] = NoiseLevel(I3);
%         sigma1 = sqrt( abs(PCANoiseLevelEstimator(I1) ));
%         sigma2 = sqrt( abs(PCANoiseLevelEstimator(I2) ));
%         sigma3 = sqrt( abs(PCANoiseLevelEstimator(I3) ));
        disp([ImgName num2str(sigma1) ' ' num2str(sigma3) ' ' num2str(sigma3)])
        fprintf(fp, '%s %f %f %f \n', [imgDir(i).name], sigma1, sigma2, sigma3);
    end
end

fclose(fp);