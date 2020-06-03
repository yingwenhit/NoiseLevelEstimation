function NLE = NLPM_NLE(percent,N_LSL,N_count)
%%% N_LSL, N_count :Local Smooth Level thresholding paramaters
 %%% precent: Gray Level thresholding paramaters
 
load('In');                  %%%%
load('Im');                  %%%%
[Nx,Ny] = size(In);         %%%%%

%%% Step2: Diffused Noise Analysis
DN = In - Im; nhood = [3 3];
localMean = filter2(ones(nhood),DN ) / prod(nhood);
localVar = filter2(ones(nhood), DN.^2) / prod(nhood) - localMean.^2;
localStd = sqrt(localVar);
figure,imshow(localStd,[]);


%%% Step3: Sepration Upon Local Smooth Level
N_GL=255;

Sigma = zeros(N_count,N_GL);
localStd1 = uni(localStd)*N_LSL; %%% modulate localStd to [0,N_LSL]
Im1 = uni(Im)*N_GL; %%% modulate Im to [0,N_GL]
for i=1:N_count
    Cor_LSL = (localStd1>N_LSL-i & ...
               localStd1<=N_LSL+1-i)*1 ...
              + zeros(Nx,Ny);
    Cor_LSL_sum = sum(Cor_LSL(:));
    Sep_LSL = Cor_LSL .* In;
    Pre_GL =uni(Cor_LSL .*Im).*Cor_LSL*N_GL;
% figure,imshow(Sep_LSL,[]),title(['Local Smooth Level=',num2str(i)]);

%%% Step4: Sepration Upon Gray Level
for j = 1:N_GL
    Cor_GL = (Pre_GL>j-1 & Pre_GL<=j)*1 + zeros(Nx,Ny);
    Cor_GL_sum = sum(Cor_GL(:));
    if Cor_GL_sum >=percent*Cor_LSL_sum
        Sep_GL = Cor_GL.*In;
        %%%%% Compute sample variance: DX = E(X^2) -(EX)^2
        E_GL = sum(Sep_GL(:))/Cor_GL_sum;   %% EX
        Sep_GL2 = Sep_GL.^2;                %% X^2
        E_GL2 =  sum(Sep_GL2(:))/Cor_GL_sum;%% E(X^2)
        D_GL = E_GL2-(E_GL)^2;              %% DX = E(X^2) -(EX)^2
        Sigma(i,j) = sqrt(D_GL);            %% Sigma is Std
             
%         figure(256),imshow(Sep_GL,[]);
%         title(['Local Smooth Level=',num2str(i),', Grey Level=',num2str(j)]);
    end
end
end
Sigma_non0 = Sigma(Sigma~=0);
NLE = mean(Sigma_non0);



