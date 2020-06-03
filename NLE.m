% Input: f is the single noisy image
% Output? NLE is the noise level estimation (standard deviation)
% All the parameters are chosen by experiments, and fixed for actual use.

%%%%%%%%%%% There is bug, 2020/1/7 %%%%%%%%%%%%%%%%%%%

function NLE = NLE(f)
% Parameters
iter = 30; K = 3; sigmaG = 1; tau = 0.1; % Parameters in NLPM
alpha = 0.8; beta = .02; % Parameters in NLE algorithm (see paper)
% NLPM preprocessing
I=f;
for i=1:iter
    I_N=I([1 1:end-1],:)-I;
    I_S=I([2:end end],:)-I;
    I_E=I(:,[2:end end])-I;
    I_W=I(:,[1 1:end-1])-I;
    
    IG = imgaussian(I,sigmaG,6*sigmaG);
    IG_N=IG([1 1:end-1],:)-IG;
    IG_S=IG([2:end end],:)-IG;
    IG_E=IG(:,[2:end end])-IG;
    IG_W=IG(:,[1 1:end-1])-IG;
    
    C_N=1./(1+abs(I_N).*abs(IG_N)/K^2);
    C_S=1./(1+abs(I_S).*abs(IG_S)/K^2);
    C_E=1./(1+abs(I_E).*abs(IG_E)/K^2);
    C_W=1./(1+abs(I_W).*abs(IG_W)/K^2);
    
    I=I+tau*(C_N.*I_N+C_S.*I_S+C_E.*I_E+C_W.*I_W);
end
figure,imshow([f,I],[])
u_pde = I;
% Remove edge to obtain flat area
edge = C_N+C_S+C_E+C_W;
flat = (edge>=alpha*max(edge(:))).*u_pde;
F = (edge>=alpha*max(edge(:)));
figure,imshow(F,[])
title('flat area')
% Compute NLE on monochrome areas
flat_int = uint8(flat);
ind = 1;
for i = 1:255
    M = (flat_int==i);
    M_num = nnz(M);
    if M_num >= beta*numel(f)
        fi = M.*f;
%         figure,imshow(fi,[]),title(num2str(i))
        fi_sum = sum(fi(:));
        fi_m = fi_sum/M_num;
        fi_s_temp = M.*(fi-fi_m).^2;
        fi_s(ind) = sum(fi_s_temp(:))/(M_num-1);
        ind = ind+1;
    end    
end

NLE = sqrt(mean(fi_s));
