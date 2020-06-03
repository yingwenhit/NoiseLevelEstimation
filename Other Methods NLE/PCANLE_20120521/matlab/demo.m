%%%%%%%%%%%%%%%%% experiments test, 2020/1/7 %%%%%%%%%%%%%%%%%%%%%%%

imgDir  = dir('../../../NLE_image/*-In-*.mat');

name_f = 'PCANLE2012_results.txt';
fp = fopen(name_f,'w');

for i = 1:length(imgDir)
    
    load(['../../../NLE_image/' imgDir(i).name]);
    x = In;
    
    if size(x,3) ~= 1
        error( 'Only grayscale 2D images can be processed.' );
    end
    
    sigma = sqrt( PCANoiseLevelEstimator(x) );
    fprintf(fp, '%f, %s \n', sigma, [imgDir(i).name]);
end

fclose(fp);