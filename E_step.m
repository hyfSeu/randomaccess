function [u,sigma,alpha] = E_step(H,y,iterNum,alpha0,a,b,delta)
%% delta为初始值取0的概率
    [M,Nc] = size(H);
    L = size(y,2);
    u = zeros(Nc,L);
    a0 = 1;
    b0 = 100;
    alpha = b0*ones(Nc,1);
    sigma = eye(Nc);
    sd = 75;          %sd为高斯分布与gamma分布的分界点，是均匀分布的下界
    su = 1e10;        %su为均匀分布上界
    eta = 1e-2;       %eta为以1为均值高斯分布的方差
    %y = y(:);
    for i = 1:iterNum
        A = diag(alpha);
        sigma = inv(alpha0*(H'*H)+A);
        u = alpha0*sigma*H'*y; %alpha0*(kron(sigma,eye(L)))*(kron(H,eye(L)))'*y;
        
        %% 2021.12.09修改版本，给超参数α赋以稀疏的先验分布
        %% 2021.12.10修改版本，α的先验从“高斯+均匀分布”改为“高斯+高斯”分布，其中另一高斯分布的均值为[100,1e10]的均匀分布
        a1 = -2/eta;
        b1 = 2/eta-diag(sigma)+diag(u*u')/L;
        c1 = 1;
        t1 = (-b1-sqrt(abs(b1).^2-4*a1*c1))/(2*a1);    %二次方程求根公式
        t2 = (1+2*a)./(diag(sigma)+diag(u*u')/L+2*b);
        t2(t2<sd) = sd;
        t2(t2>su) = su;
        t1(t1>sd) = sd;
        Q1 = 0.5*log(t1)-0.5*(diag(sigma)+diag(u*u')/L).*t1+log((1-delta));
        Q2 = 0.5*log(t2)-0.5*(diag(sigma)+diag(u*u')/L).*t2+log(delta);
        t1(Q1<Q2) = t2(Q1<Q2);
        alpha = t1;
%         alpha = (L+L*a)./(L*diag(sigma)+diag(u*u')+L*b);
%         %%初始版本，假设超参数α服从gamma分布
    end
end