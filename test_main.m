clear;
clc;
rng(175);

%% 系统初始化
K = 16;                                         %总用户数
M = 300;                                         %天线数
delta = 0.25;
% Ka = floor(0.125*K*M);                                         %活跃用户数
L = 6;                                          %信号长度
Lk = 25;                                         %信号路径数量最大值
pl_exp = 3.7;                                       %路径损耗因子
d0 = 10;
pmax = 17;
MC_num = 4;                                   %monte-calo 仿真次数
current_error = zeros(MC_num,1);
% user_index = randperm(K*M,Ka);                     %随机获取用户索引
% user_index = sort(user_index);
% theta = 2*pi*rand(Nc,1);
estimated_error = [];
error = [];
cur_date = date;
cur_time = fix(clock);
str = sprintf('%s %.2d:%.2d:%.2d\n',cur_date,cur_time(4),cur_time(5),cur_time(6));
disp(str);
for i = 1:MC_num
    distance = 250*rand(K,1);                       %各用户到基站的距离，方圆1Km
    LS = get_snr_dB(distance,d0,pmax,pl_exp);
    lamda = 10.^(LS/10);                 %(distance.^mu).*10^13.4;                 %大尺度参数，这个在基站是已知的
    % lamda =1e4*ones(K,1);
    % G = (randn(K,M)+1j*randn(K,M))/sqrt(2);     %实际的信道矩阵
    % G = G';

    supp = zeros(K,M);                               %支持矩阵，表面接收信道矩阵非零元素所在，1表示活跃，0表示静默
    A = gen_Unit(M);                                 %生成调制酉矩阵
    % supp(user_index,:) = ones(Ka,M);
    % tempidx = ceil(Lk*rand(Ka,1))+3; 
    % for i = 1:Ka
    %     rowidx = floor(user_index(i)/M)+1;
    %     colidx = mod(user_index(i),M);
    %     if(colidx==0)
    %         colidx = M;
    %     end
    %     supp(rowidx,colidx)=1;  %支持向量
    % end
    % G = G.*supp(1:K,1:M);
    message = (randsrc(K,M)+1)/2;
    G = sp_coding(message,delta);
    [~,tempM] = size(G);
    [trueK,~] = size(unique(G,'rows'));
    %G(:,1) = (ones(K,1)+1j*ones(K,1))/sqrt(2);
    % delta = sum(supp(:))/(K*M);                      %信道稀疏度
    X = (randn(L,K)+1j*randn(L,K))/sqrt(2);
    X(1,:) = (ones(1,K)+1j*ones(1,K))/sqrt(2);
    X = X*sqrt(diag(lamda));
    signals = X*G;                                   %接收信号为信号乘以信道
    noise1 = wgn(L,tempM,0,'complex');                   %噪声信号，功率归一化

    %% 接收信号
    Y = signals + noise1;
    noise_power = sum(sum(abs(noise1).^2))/(M*L);
    signal_power = sum(signals'*signals)/(M*L);
    %Z = Y*A*inv(A'*A);

    %% 循环迭代
    % initialize Xhat,Hhat,Xvar,Hvar,Shat
    Lamda = repmat(lamda',L,1);                        %var of channels
    Xmu = Lamda;
    Xmean = zeros(L,K);
    Hmu = 0.01*ones(K,tempM);
    Hmean = ones(K,tempM);
    
    [Xhat,Xvar,Hhat,Hvar] = BiG_AMP(Y,noise_power,K,delta,Xmean,Xmu,Hmean,Hmu,[],[],[],[]);
    [orignal_code,hamming_distance,corres_idx] = sp_decoding(Hhat,message,delta);
    idx1 = unique(corres_idx);
    current_error(i) = (sum(hamming_distance*M)+M*(trueK-length(idx1)))/(trueK*M);
end
estimated_error = [estimated_error,current_error];
error = mean(current_error);
cur_date = date;
cur_time = fix(clock);
str = sprintf('%s %.2d:%.2d:%.2d\n',cur_date,cur_time(4),cur_time(5),cur_time(6));
disp(str);
% num_Message = 100:50:1000;
% plot(num_Message,error,'k');hold on
% plot(num_Message,error,'ks');
% xlabel('length of messages','FontSize',10.508);
% ylabel('error_rate','FontSize',10.508)