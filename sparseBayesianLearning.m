clear;
clc;
Nc = 512;          %总用户数
Ka = 13;            %活跃用户数
M = 120;             %天线数
L = 128;            %信号长度
mu = -3.7;          %路径损耗因子
a = 1e-8;
b = 1e-8;
c = 1e-8;
d = 1e-8;
user_index = randperm(Nc,Ka);    
user_index = sort(user_index);
theta = 2*pi*rand(Nc,1);
distance = 1000*rand(Nc,1);      %各用户到基站的距离，方圆1Km
%lamda = 150^mu*ones(Nc,1);
% snr = sys_snr(17,10);
% lamda = (2*10^-5.45)./(1+(distance/10).^(-mu));            %大尺度系数
lamda1 = distance.^mu;
sk = 10*log10(lamda1)+121;
%lamda = 17+104+10*log10(distance.^mu);
lamda = distance.^mu;
G = sqrt(diag(lamda))*(randn(Nc,M)+1j*randn(Nc,M))/sqrt(2);     %实际的信道矩阵
G = G';
% G = 1e-3*(randn(M,Nc)+1j*randn(M,Nc))/sqrt(2);
message = zeros(Nc,L);           %活跃用户数据消息预设为随机数值，其余为0
message(user_index,:)=((1-2*fix(rand(Ka,L)./0.5))+1j*(1-2*fix(rand(Ka,L)./0.5)))/sqrt(2);
signals = G*message;
%y = awgn(signals,15,'measured');
noise1 = wgn(M,L,-134,'complex');
%noise2 = y-signals;
y = signals + noise1;
noise_power = sum(sum(abs(noise1).^2))/(M*L);
signal_power = sum(signals'*signals)/(M*L);
%y = signals+noise1;
%db = 10*log10(signal_power/noise_power);

%% 先假设H=G，用以验证稀疏贝叶斯学习的有效性
H = G;
iterNum = 500;
alpha = ones(Nc,L);
u = zeros(Nc,L);
sigma = zeros(Nc,Nc,L);
%alpha0 = 1;
%alpha = 5e-11*(randn(Nc,L)+1j*randn(Nc,L))/sqrt(2);
for i = 1:L
    alpha0 = 1/noise_power;
    %alpha0 = 1;
    A = diag(alpha(:,i));
    sigma(:,:,i) = inv(alpha0*(H'*H)+A);
    u(:,i) = alpha0*sigma(:,:,i)*H'*y(:,i);
    for j =1:iterNum
        gamma = 1-alpha(:,i).*diag(sigma(:,:,i));
        %alpha(:,i) = gamma./(u.*conj(u));
        alpha(:,i) = (1+2*a)./((u(:,i).*conj(u(:,i)))+diag(sigma(:,:,i))+2*b);
        %alpha0 = (M-sum(gamma))/((y(:,i)-H*u)'*(y(:,i)-H*u));
        %alpha0 = (M+2*c)/((y(:,i)-H*u)'*(y(:,i)-H*u)+sum(gamma)/alpha0+2*d);
        A = diag(alpha(:,i));
        sigma(:,:,i) = inv(alpha0*(H'*H)+A);
        u(:,i) = alpha0*sigma(:,:,i)*H'*y(:,i);
    end
end


tt=(y*u')*inv(sum(sigma,3)+u*u');
temp1=u(user_index);
temp2=message(user_index);
