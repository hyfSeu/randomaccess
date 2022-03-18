clear;
clc;

%% 系统初始化
K = 256;                                         %总用户数
M = 64;                                         %天线数
Ka = floor(0.1*K*M);                                         %活跃用户数
L = 128;                                          %信号长度
Lk = 25;                                         %信号路径数量最大值
mu = -3.7;                                       %路径损耗因子
beta = 0.05;                                        %damping系数
tao = 1e-6;                                      %迭代终止阈值
IterNum = 1500;                                  %AMP算法迭代次数
IterMin = 30;
stepWindow = 1;                                 %size of stepWindw
valOpt = [];                                     %set of succeed step results
stepTol = -1;
stepMin = 0.05;
stepMax = 0.5;
stepIncr = 1.1;
stepDecr = 0.5;
maxStepDecr = 0.5;   
% delta = (K-Ka)/K;                              %用户稀疏程度
% a = 1e-10;                                     %a,b,c,d为gamma分布参数
% b = 1e-10;
% c = 1e-8;
% d = 1e-8;
PvarMin = 1e-13;                                    %% warning!!!
ZvarToPvarMax = 0.99;                              %% warning!!!
varThresh = 1e6;                                   %% warning!!!
HvarMin = 0;                                    %% warning!!!
XvarMin = 0;                                    %% warning!!!
user_index = randperm(K*M,Ka);                     %随机获取用户索引
user_index = sort(user_index);
% theta = 2*pi*rand(Nc,1);
distance = 1000*rand(K,1);                       %各用户到基站的距离，方圆1Km
lamda = 100*ones(K,1);                 %(distance.^mu).*10^13.4;                 %大尺度参数，这个在基站是已知的
G = sqrt(diag(lamda))*(randn(K,M)+1j*randn(K,M))/sqrt(2);     %实际的信道矩阵
% G = G';
supp = zeros(K,M);                               %支持矩阵，表面接收信道矩阵非零元素所在，1表示活跃，0表示静默
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
signals = X*G;                                   %接收信号为信号乘以信道
noise1 = wgn(L,M,-20,'complex');                   %噪声信号，功率归一化




%% 接收信号
Y = signals + noise1;
noise_power = sum(sum(abs(noise1).^2))/(M*L);
signal_power = sum(signals'*signals)/(M*L);




%% 循环迭代
% initialize Xhat,Hhat,Xvar,Hvar,Shat
Lamda = repmat(lamda,1,M);                        %var of channels
Pvar_mean = 1;
valIn = 0;
val = zeros(IterNum,1);
failCount = 0;
step = 1;

%X variables
Xhat = (randn(L,K)+1j*randn(L,K))/sqrt(2);
Xvar = 10*ones(L,K);
XhatBar = zeros(L,K);
XhatOpt = zeros(L,K);
XhatBarOpt = zeros(L,K);

%H variables
Hhat = zeros(K,M);
Hvar = 10*Lamda;
HhatBar = zeros(K,M);
HhatOpt = zeros(K,M);
HhatBarOpt = zeros(K,M);

%S variables
shat = zeros(L,M);
svar = zeros(L,M);
shatOpt = zeros(L,M);
svarOpt = zeros(L,M);
shatNew = zeros(L,M);
svarNew = zeros(L,M);

%Variance momentum terms
PvarOpt = zeros(L,M);
ZvarOpt = zeros(L,M);
ZhatOpt = zeros(L,M);

it = 0;
stop = false;

%main itereration loop
while ~stop
   it = it+1;
   if(it >= IterNum)
       stop = true;
   end
   if(it == 1)
       step = 1;
   end
   Xhat2 = abs(Xhat).^2;
   Hhat2 = abs(Hhat).^2;
   Zvar = Xvar*Hhat2+Xhat2*Hvar;
   Pvar = Zvar+Xvar*Hvar;
   Zvar = step*Zvar+(1-step)*ZvarOpt;
   Pvar = step*Pvar+(1-step)*PvarOpt;
   Zhat = Xhat*Hhat;
   Phat = Zhat-shat.*(Zvar/Pvar_mean);
   
   %compute log likelihood at the output and add it the total negative
   %K-L distance at the input.
   valOut = sum(sum(logLike(Y,noise_power,Zhat,Pvar)));
   val(it) = valOut + valIn;
   
   %determine if candidate passed
   if ~isempty(valOpt)
       stopInd = length(valOpt);
       startInd = max(1,stopInd - stepWindow);
       pass = (val(it)>min(valOpt(startInd:stopInd)))||(beta<=stepMin);
   else
       pass = true;
   end
   
   % If pass, set the optimal values and compute a new target shat and
   % snew.
   if pass
       beta = stepIncr*beta;
       
       shatOpt = shat;
       svarOpt = svar;
       XhatBarOpt = XhatBar;
       XhatOpt = Xhat;
       HhatBarOpt = HhatBar;
       HhatOpt = Hhat;
       PvarOpt = Pvar;
       ZvarOpt = Zvar;
       
       %Bound pvar
       Pvar = max(Pvar,PvarMin);
       
       %only keep the successful step valOpt values
       valOpt = [valOpt val(it)];
       
       %compute the nonlinear step
       Zhat0 = (noise_power*Phat+Y.*Pvar)./(noise_power+Pvar);
       Zvar0 = (noise_power*Pvar)./(noise_power+Pvar);
       
       %compute 1/pvar
       PvarInv = Pvar_mean./Pvar;
       
       %update the shat quantities
       shatNew = PvarInv.*(Zhat0-Phat);
       svarNew = PvarInv.*(1-min(Zvar0./Pvar,ZvarToPvarMax));
       
       %enforce step size bounds
       beta = min([max([beta stepMin]) stepMax]);
       
%        XvarOpt = Xvar;
%        HvarOpt = Hvar;
%        ZhatOpt = Zhat;
%    
%        %caclulate the h variables
%        Kappa = 1./(1+(1-delta)/delta*sqrt((Lamda+Rvar)./Rvar).*exp(-Lamda.*(abs(Rhat).^2)./(2*Rvar.*(Rvar+Lamda))));
%        Hhat = Kappa.*(Lamda.*Rhat)./(Lamda+Rvar);
%        Hvar = Kappa.*(abs((Lamda.*Rhat)./(Lamda+Rvar)).^2+Lamda.*Rvar./(Lamda+Rvar))-abs(Hhat).^2;
% 
%        %caclulate the x variables
%        Xhat = Qhat./(1+Qvar);                       %px is a Normal distrbution
%        Xvar = Qvar./(1+Qvar);
%        beta = min([max([beta stepMin]) stepMax]);
   else
       beta = stepDecr*beta;
       beta = max(stepMin,beta);
       if(beta<stepTol)
           stop = true;
       end
   end
   
   %check for convergence if step was successful
   if pass
       if any(isnan(Zhat(:)))||any(isinf(Zhat(:)))
           stop = true;
       else
           testVal = norm(Zhat(:)-ZhatOpt(:))/norm(Zhat(:));
           if (it>1)&&(testVal<tao)
               stop = true;
           end
       end
       XvarOpt = Xvar;
       HvarOpt = Hvar;
       ZhatOpt = Zhat;
   end
   
   if(it>1)
       step = beta;
   end
   shat = (1-step)*shatOpt+step*shatNew;
   svar = (1-step)*svarOpt+step*svarNew;
   HhatBar = (1-step)*HhatBarOpt+step*HhatOpt;
   XhatBar = (1-step)*XhatBarOpt+step*XhatOpt;
   
   %caclulate the r variables
   Rvar = 1./((abs(XhatBar).^2)'*svar);
   Rvar(Rvar>varThresh) = varThresh;
   RGain = (1-(Rvar.*(Xvar'*svar)));
   RGain = min(1,max(0,RGain));
   Rhat = HhatBar.*RGain+Rvar.*(XhatBar'*shat);
   Rvar = max(Rvar,HvarMin);
   
   %caclulate the q variables
   Qvar = 1./(svar*(abs(HhatBar).^2)');
   Qvar(Qvar>varThresh) = varThresh;
   QGain = (1-(Qvar.*(svar*Hvar')));
   QGain = min(1,max(0,QGain));
   Qhat = XhatBar.*QGain+Qvar.*(shat*HhatBar');
   Qvar = max(Qvar,XvarMin);
   
   %caclulate the h variables
   Kappa = 1./(1+(1-delta)/delta*sqrt((Lamda+Rvar)./Rvar).*exp(-Lamda.*(abs(Rhat).^2)./(2*Rvar.*(Rvar+Lamda))));
   Hhat = Kappa.*(Lamda.*Rhat)./(Lamda+Rvar);
   Hvar = Kappa.*(abs((Lamda.*Rhat)./(Lamda+Rvar)).^2+Lamda.*Rvar./(Lamda+Rvar))-abs(Hhat).^2;
   temp = 0.5*(log(Rvar./(Rvar+Lamda))+Lamda./(Rvar+Lamda)-abs(Lamda.*Rhat./(Rvar+Lamda)).^2./Lamda);
   valInH = Kappa.*temp+Kappa.*log(max(1e-8,delta)./max(Kappa,1e-8))+(1-Kappa).*log(max(1e-8,(1-delta))./max((1-Kappa),1e-8));

   %caclulate the x variables
   Xhat = Qhat./(1+Qvar);                       %px is a Normal distrbution
   Xvar = Qvar./(1+Qvar);
   valInX = log(Xvar)+(1-Xvar)-abs(Xhat).^2;    %specific
    
   %update valIn
   valIn = sum(valInH(:))+sum(valInX(:));
   
   %don't stop before minimum iteration count
   if it<IterMin
       stop = false;
   end
       
   
%    Pvar = max(Pvar,PvarMin);
%    PvarInv = 1./Pvar;
%    Zhat0 = (noise_power*Phat+Y.*Pvar)./(noise_power+Pvar);
%    Zvar0 = (noise_power*Pvar)./(noise_power+Pvar);
%    shatNew = PvarInv.*(Zhat0-Phat);
% %    svarNew = PvarInv.*(1-Zvar0./Pvar);
%    svarNew = PvarInv.*(1-min(Zvar0./Pvar,ZvarToPvarMax));
   
%    testVal = norm(Zhat(:)-ZhatOpt(:))/norm(Zhat(:));
%    
%    
%    if(beta<stepTol)
%        break;
%    end
%    
%    shat = (1-beta)*shatOpt+beta*shatNew;
%    svar = (1-beta)*svarOpt+beta*svarNew;
%    HhatBar = (1-beta)*HhatBarOpt+beta*HhatOpt;
%    XhatBar = (1-beta)*XhatBarOpt+beta*XhatOpt;
   

   
   



%    if pass
%        beta = beta+stepIncr;
%        valOpt = [valOpt costJ];
%        XvarOpt = Xvar;
%        HvarOpt = Hvar;
%        ZhatOpt = Zhat;
%    
%        %caclulate the h variables
%        Kappa = 1./(1+(1-delta)/delta*sqrt((Lamda+Rvar)./Rvar).*exp(-Lamda.*(abs(Rhat).^2)./(2*Rvar.*(Rvar+Lamda))));
%        Hhat = Kappa.*(Lamda.*Rhat)./(Lamda+Rvar);
%        Hvar = Kappa.*(abs((Lamda.*Rhat)./(Lamda+Rvar)).^2+Lamda.*Rvar./(Lamda+Rvar))-abs(Hhat).^2;
% 
%        %caclulate the x variables
%        Xhat = Qhat./(1+Qvar);                       %px is a Normal distrbution
%        Xvar = Qvar./(1+Qvar);
%        beta = min([max([beta stepMin]) stepMax]);
%    else
%        beta = beta-stepDecr;
%        beta = max(stepMin,beta);
%        continue;
%    end
end
