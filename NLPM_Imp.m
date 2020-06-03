% Apply AOS to solve semi-implicit scheme
function u = NLPM_Imp(f, iter, K, dt, sigma)

[Nx, Ny] = size(f);
% dt = 10;
% K = 8;
K = K^2;
% iter = 1;
% parameter for gaussian blur: sigma, siz
% sigma = 1;
siz = 6*sigma;
fg = imgaussian(f, sigma, siz);

[Gfgx, Gfgy] = gradient(fg);
Gfg = sqrt(Gfgx.^2 + Gfgy.^2);

u = f;
%% AOS Algorithm
for i=1:iter
    [Gux, Guy] = gradient(u);
    Gu = sqrt(Gux.^2 + Guy.^2);
    g = 1./(1+Gfg.*Gu/K);
    D_c = u;
    
%% Column
    
    Alpha_c = 1 + dt*(g([1, 1:end-1],:) + 2*g + g([2:end end],:))/2;
    Alpha_c(1,:) = 1 + dt*(g(1,:) + g(2,:))/2;
    Alpha_c(end,:) = 1 + dt*(g(end,:) + g(end-1,:))/2;
    G = -1*dt*(g([1:end end],:) + g([1, 1:end],:))/2;
    Beta_c = G([2:end-1],:);
    Gamma_c = G([2:end-1],:);
%     u = Thomas(Alpha_c, Beta_c, Gamma_c, D_c);
    for j=1:Nx
        Alpha = Alpha_c(:,j);
        Beta  = Beta_c(:,j);
        Gamma = Gamma_c(:,j);
        D     = D_c(:,j);
        u_c(:,j) = Thomas(Alpha, Beta, Gamma, D);
        % Thomas for row
        % 
    end
%% Row
    Alpha_r = 1 + dt*(g(:,[1, 1:end-1]) + 2*g + g(:,[2:end end]))/2;
    Alpha_r(:,1) = 1 + dt*(g(:,1) + g(:,2))/2;
    Alpha_r(:,end) = 1 + dt*(g(:,end) + g(:,end-1))/2;
    G = -1*dt*(g(:,[1:end end]) + g(:,[1, 1:end]))/2;
    Beta_r = G(:,[2:end-1]);
    Gamma_r = G(:,[2:end-1]);
%     u = Thomas(Alpha_r, Beta_r, Gamma_r, D_c);
    for j=1:Ny
        Alpha = Alpha_r(j,:);
        Beta  = Beta_r(j,:);
        Gamma = Gamma_r(j,:);
        D     = D_c(j,:);
        u_rr(:,j) = Thomas(Alpha, Beta, Gamma, D);
        % Thomas for column
    end
    u_r = u_rr';
    
    u = (u_c + u_r)/2;
%     figure(111),imshow([u,f],[]);title("Thomas Implicit");

end

end

function X = Thomas(Alpha, Beta, Gamma, D)

%% Check

N  = numel(Alpha);
N_  = numel(Beta);
N__ = numel(Gamma);
N___= numel(D);
if ((N_~=N__) && (N~=N_-1) && (N~=N___))
    error('Dimension error!');
end
M = zeros(N,1);
L = zeros(N-1,1);
X = zeros(N,1);


%% RL decomposition
R = Beta;
M(1) = Alpha(1);
for i=1:N-1
    L(i) = Gamma(i)./M(i);
    M(i+1) = Alpha(i+1) - L(i).*Beta(i);
end

%% Forward
Y = zeros(N,1);
Y(1) = D(1);
for i=2:N
    Y(i) = D(i) - L(i-1).*Y(i-1);
end

%% Backward
X(N) = Y(N)/M(N);
for i=N-1:-1:1
    X(i) = (Y(i) - Beta(i).*X(i+1))/M(i);
end

end