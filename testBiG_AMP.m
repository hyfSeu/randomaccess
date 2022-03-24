clear;
clc;
rng(676);

%% 系统初始化
K = 128;                                         %总用户数
M = 3200;                                         %天线数
Ka = floor(0.2*K*M);                                         %活跃用户数
L = 64;                                          %信号长度
Lk = 25;                                         %信号路径数量最大值
mu = -3.7;                                       %路径损耗因子
user_index = randperm(K*M,Ka);                     %随机获取用户索引
user_index = sort(user_index);
% theta = 2*pi*rand(Nc,1);
distance = rand(K,1);                       %各用户到基站的距离，方圆1Km
lamda = ones(K,1);                 %(distance.^mu).*10^13.4;                 %大尺度参数，这个在基站是已知的
G = sqrt(diag(lamda))*(randn(K,M)+1j*randn(K,M))/sqrt(2);     %实际的信道矩阵
% G = G';
supp = zeros(K,M);                               %支持矩阵，表面接收信道矩阵非零元素所在，1表示活跃，0表示静默
% supp(user_index,:) = ones(Ka,M);
% tempidx = ceil(Lk*rand(Ka,1))+3; 
for i = 1:Ka
    rowidx = floor(user_index(i)/M)+1;
    colidx = mod(user_index(i),M);
    if(colidx==0)
        colidx = M;
    end
    supp(rowidx,colidx)=1;  %支持向量
end
G = G.*supp;
delta = sum(supp(:))/(K*M);                      %信道稀疏度
X = (randn(L,K)+1j*randn(L,K))/sqrt(2);
X(1,:) = (ones(1,K)+1j*ones(1,K))/sqrt(2);
signals = X*G;                                   %接收信号为信号乘以信道
noise1 = wgn(L,M,-20,'complex');                   %噪声信号，功率归一化

%% 接收信号
Y = signals + noise1;
noise_power = sum(sum(abs(noise1).^2))/(M*L);
signal_power = sum(signals'*signals)/(M*L);

%% 循环迭代
% initialize Xhat,Hhat,Xvar,Hvar,Shat
Lamda = repmat(lamda,1,M);                        %var of channels
Xmu = 1;
Hmu = Lamda;
[Xhat,Xvar,Hhat,Hvar] = BiG_AMP(Y,noise_power,K,delta,Xmu,Hmu);
