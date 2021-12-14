clear;
clc;
Nc = 512;          %总用户数
Ka = 13;            %活跃用户数
M = 120;             %天线数
L = 32;            %信号长度
mu = -3.7;          %路径损耗因子
EMstep = 300;       %EM算法迭代次数
delta = (Nc-Ka)/Nc; 
a = 1e-10;           %a,b,c,d为gamma分布参数
b = 1e-10;
c = 1e-8;
d = 1e-8;
user_index = randperm(Nc,Ka);      %随机获取用户索引
user_index = sort(user_index);
theta = 2*pi*rand(Nc,1);
distance = 1000*rand(Nc,1);      %各用户到基站的距离，方圆1Km
lamda = (distance.^mu).*10^13.4;                   %大尺度参数
G = sqrt(diag(lamda))*(randn(Nc,M)+1j*randn(Nc,M))/sqrt(2);     %实际的信道矩阵
G = G';
% G = 1e-3*(randn(M,Nc)+1j*randn(M,Nc))/sqrt(2);
message = zeros(Nc,L);           %活跃用户数据消息预设为随机数值，其余为0
% message(user_index,:)= (randn(Ka,L)+1j*randn(Ka,L))/sqrt(2);
message(user_index,:)= [(ones(Ka,1)+1j*ones(Ka,1))/sqrt(2),((1-2*fix(rand(Ka,L-1)./0.5))+1j*(1-2*fix(rand(Ka,L-1)./0.5)))/sqrt(2)];
% message(user_index,:) = repmat( [(ones(Ka,1)+1j*ones(Ka,1))/sqrt(2),(ones(Ka,1)-1j*ones(Ka,1))/sqrt(2),(-ones(Ka,1)-1j*ones(Ka,1))/sqrt(2),(-ones(Ka,1)+1j*ones(Ka,1))/sqrt(2)],1,L/4);
signals = G*message;
%y = awgn(signals,15,'measured');
noise1 = wgn(M,L,0,'complex');
%noise2 = y-signals;
y = signals + noise1;
noise_power = sum(sum(abs(noise1).^2))/(M*L);
signal_power = sum(signals'*signals)/(M*L);
%y = signals+noise1;
%db = 10*log10(signal_power/noise_power);
iterNum = 500;
alpha0 = 1/noise_power;

% [u,sigma,alpha] = E_step(G,y,iterNum,alpha0,a,b);
% H = (y*u')*inv(sum(sigma,3)+u*u');
% t = message(user_index,:);
% t0 = x(user_index,:);

%% 信道矩阵初始化
distance_rand = 1000*rand(Nc,1);
lamda_rand = (distance_rand.^mu).*10^13.4;
H = sqrt(diag(lamda_rand))*(randn(Nc,M)+1j*randn(Nc,M))/sqrt(2);
H = H';

%% 以下代码用以测试稀疏信号检测算法
% [r1,u_r1,x,alpha1] = VAMP(G,y,iterNum,1,delta,noise_power);
%[r,u_r,x,sigma] = New_AMP(G,y,iterNum,1,delta,noise_power);
% [u,sigma,alpha] = E_step(G,y,iterNum,alpha0,a,b,delta);
% t1 = u(user_index,:);
% t2 = message(user_index,:);
% H = (y*u')*inv(L*sigma+u*u');

% x = 1e-8*ones(Nc,L);
% user_index0 = randperm(Nc,Ka);      %随机获取用户索引
% user_index0 = sort(user_index0);
% x(user_index0,:) = repmat( [(ones(Ka,1)+1j*ones(Ka,1))/sqrt(2),(ones(Ka,1)-1j*ones(Ka,1))/sqrt(2),(-ones(Ka,1)-1j*ones(Ka,1))/sqrt(2),(-ones(Ka,1)+1j*ones(Ka,1))/sqrt(2)],1,L/4);
% sigma = 1e-8*eye(Nc,Nc);
% H0 = (y*x')*inv(L*sigma+diag(x*x'));
% H(:,user_index0) = H0(:,user_index0);
% [u,sigma,alpha] = E_step(H,y,iterNum,alpha0,a,b,delta);

%% EM算法过程
% for i = 1:EMstep
%     [r,u_r,x,sigma] = AMP(H,y,iterNum,1,delta,noise_power);
%     H = (y*x')*inv(reshape(sum(sigma,1),[Nc,Nc])+x*x');
% end
% [r,u_r,x,sigma] = AMP(G,y,iterNum,1,delta,noise_power);
% H = (y*x')*inv(reshape(sum(sigma,1),[Nc,Nc])+x*x');
% ttt = G(:,user_index);
% ttt0 = H(:,user_index);
% [r0,u_r0,x0,u_x0] = AMP(G,y,2*iterNum,1,delta,noise_power);
% tt0 = x0(user_index,:);
% [u,sigma,alpha] = E_step(G,y,iterNum,alpha0,a,b,delta);
% ttt = u(user_index,:);
% [u1,sigma1,alpha1] = E_step(G,y,2*iterNum,alpha0,a,b);
% ttt1 = u1(user_index,:);


    
%% EM算法过程
cur_date = date;
cur_time = fix(clock);
str = sprintf('%s %.2d:%.2d:%2d\n', cur_date,cur_time(4),cur_time(5),cur_time(6));
disp(str);
disp("start!");
for i = 1:EMstep
    [u,sigma,alpha] = E_step(H,y,iterNum,alpha0,a,b,delta);
    H = (y*u')*inv(sum(sigma,3)+u*u');
    if(mod(i,EMstep/10)==0)
        disp(i/EMstep);
        disp('\n');
    end
end


% tt=(y*u')*inv(sum(sigma,3)+u*u');
% temp1=u(user_index);
% temp2=message(user_index);
