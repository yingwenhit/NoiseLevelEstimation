% 2020/05/26 Implicit Scheme for NLPM, with Gauss-Seidel Iteration
function u = NLPM_Imp_GS(f, iter, K, dt, sigmaG)

K = K^2;
% iter = 1;
% parameter for gaussian blur: sigma, siz
sigma = sigmaG;
siz = 6*sigmaG;

fg = imgaussian(f, sigma, siz);
[Gfgx, Gfgy] = gradient(fg);
Gfg = sqrt(Gfgx.^2 + Gfgy.^2);
u = f;

for i=1:iter
[Gux, Guy] = gradient(u);
Gu = sqrt(Gux.^2 + Guy.^2);
g = 1./(1+Gfg.*Gu/K);

[u] = Imp_NLPM_GS2(g, u, dt);
% figure(111),imshow([u,f],[]);title('Imp G-S')
end

end