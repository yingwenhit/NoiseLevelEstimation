% Input: f is the single noisy image
% Output: NLE is the noise level estimation (standard deviation)
% All the parameters are chosen by experiments, and fixed for actual use.
% According to the chapter 2 in paper.

% Update,2020/06/01: Using NLPM semi-implict scheme and Gauss-Seidel Iteration method
function [NLE, u_pde, O_flat, Count] = NLPM_miniBatch_Imp_NLE(f_, iter, K, sigmaG, tau, alpha, beta, Bch)


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
u_pde_ = NLPM_Imp_GS(f_, iter, K, tau, sigmaG)+1;
% u_pde = NLPM_Imp(f, iter, K, tau, sigmaG)+1;
% figure(111),imshow([u_pde,f],[]);title('Imp G-S')
[Nx, Ny] = size(u_pde_);
NLE_ = 0;
count = 1;
Ns = ceil(Nx/Bch);
Ms = ceil(Ny/Bch);
for ki=1:Ns
    for kj=1:Ms
%         disp([num2str(ki) ' ' num2str(kj)]);
        kii = (ki-1)*Bch+1;
        kjj = (kj-1)*Bch+1;
        kii_ = kii+Bch-1;
        if kii_> Nx
            kii_ = Nx;
            dif_i = Nx - kii + 1;
            if dif_i < 10
                continue;
            end
        end
        kjj_ = kjj+Bch-1;
        if kjj_>Ny
            kjj_ = Ny;
            dif_j = Ny - kjj + 1;
            if dif_j < 10
                continue;
            end
        end
%         disp([num2str(kii_) ' ' num2str(kjj_)]);
        u_pde = u_pde_([kii:kii_], [kjj:kjj_]);
        f = f_([kii:kii_],[kjj:kjj_]);
        
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
        % figure,imshow(O_flat);title('O\_flat')
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
            if (M_num >= beta_f) && (M_num>1)
                fi = M.*f;
                %         figure,imshow(fi,[]),title(num2str(i))
                fi_sum = sum(fi(:));
                fi_m = fi_sum/M_num;
                fi_s_temp = M.*(fi-fi_m).^2;

                fi_s(ind) = sum(fi_s_temp(:))/(M_num-1);
                ind = ind+1;
            end
        end
        if (exist('fi_s'))
            NLE_(count) = sqrt(mean(fi_s));
            count = count+1;
        end
    end
end
NLE = mean(NLE_);
