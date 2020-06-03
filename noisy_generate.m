clear all
close all

for std_n = [10, 15, 20, 25]
    I = imread('../Image/mdn4.bmp');
    if (ndims(I) == 3)
        I = rgb2gray(I);
    end
    % I = double(rgb2gray(I));
    I = double(I);

%     std_n = 30;
    NI = randn(size(I))*std_n; % White Gaussian noise
    In = I + NI;  % noisy input image
    figure,imshow(In,[])
    save (['mdn4_w-In-' num2str(std_n)], 'In', 'I', 'std_n')
end