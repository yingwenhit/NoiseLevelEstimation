% Input: f is the single noisy image
% Output: NLE is the noise level estimation (standard deviation)
% All the parameters are chosen by experiments, and fixed for actual use.
% According to the chapter 2 in paper.

% Update,2020/06/01: Using NLPM semi-implict scheme and Gauss-Seidel Iteration method
function [NLE, u_pde, O_flat, Count] = NLPM_Imp_NLE(f, iter, K, sigmaG, tau, alpha, beta)


% NLPM preprocessing
% I=f;
% for i=1:iter
%     I_N=I([1 1:end-1],:)-I;
%     I_S=I([2:end end],:)-I;
%     I_E=I(:,[2:end end])-I;
%     I_W=I(:,[1 1:end-1])-I;
%     
%     IG = imgaussian(I,sigmaG,6*sigmaG);
%     IG_N=IG([1 1:end-1],:)-IG;
%     IG_S=IG([2:end end],:)-IG;
%     IG_E=IG(:,[2:end end])-IG;
%     IG_W=IG(:,[1 1:end-1])-IG;
%     
%     C_N=1./(1+abs(I_N).*abs(IG_N)/K^2);
%     C_S=1./(1+abs(I_S).*abs(IG_S)/K^2);
%     C_E=1./(1+abs(I_E).*abs(IG_E)/K^2);
%     C_W=1./(1+abs(I_W).*abs(IG_W)/K^2);
%     
%     I = I + tau*(C_N.*I_N+C_S.*I_S+C_E.*I_E+C_W.*I_W);
% %     figure(33),imshow(uint8([f,I]));title(['i=' num2str(i)])
% end
% u_pde = I+1;
u_pde = NLPM_Imp_GS(f, iter, K, tau, sigmaG)+1;
% u_pde = NLPM_Imp(f, iter, K, tau, sigmaG)+1;
% figure(111),imshow([u_pde,f],[]);title('Imp G-S')

% Remove edge to obtain flat area
G_u_pde = imgaussian(u_pde,sigmaG,6*sigmaG);
[grad_Gx, grad_Gy] = gradient(G_u_pde);
[ux, uy] = gradient(u_pde);
Grad1 = sqrt(grad_Gx.^2 + grad_Gy.^2);
Grad2 = sqrt(ux.^2 + uy.^2);
Grad = Grad1.*Grad2;

G_inf = 1./(1 + Grad);

O_flat = (1 + Grad)<=(1/alpha);
flat = double(O_flat) .* u_pde;
figure,imshow(O_flat);title('O\_flat')
% edge = C_N+C_S+C_E+C_W;
% flat = (edge>=alpha*max(edge(:))).*u_pde;

% Compute NLE on monochrome areas
flat_int = uint8(flat);
ind = 1;    
for i = 1:256
    M = (flat_int == i);
    M_num = nnz(M);
    Count(i) = M_num;
end

Max_num = max(Count);
beta_f = beta*numel(f);

if beta_f >= Max_num
    beta_f = Max_num/2;
end

for i = 1:256
    M = (flat_int == i);
    M_num = Count(i);
    if M_num >= beta_f
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
