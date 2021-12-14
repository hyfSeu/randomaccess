
function [r,u_r,x,u_x] = New_AMP(H,y,iterNum,mu,lamda,mu0)
%% H是信道矩阵，y是接收信号，iterNum是迭代次数，alpha为信号高斯分布均值，
%% mu为信号高斯分布方差，lamda为伯努利分布参数，mu0为噪声信号方差
    [M,Nc] = size(H);
    L = size(y,2);
    x = zeros(Nc,L) + 0.023;
    u_x = 0.09+zeros(Nc,1);
    s = zeros(M,L);
    u_s = zeros(M,1);
    r = zeros(Nc,L);
    u_r = zeros(Nc,1);
    p = zeros(M,L);
    %sigma = zeros(L,Nc,Nc);
    for i = 1:iterNum
        u_p = (abs(H).^2)*u_x;
        p = H*x-diag(u_p)*s;
        z = H*x;
        s = diag(1./(u_p+mu0))*(y-p);
        u_s = 1./(u_p+mu0);
        u_r = 1./((abs(H).^2)'*u_s);
        r = x+diag(u_r)*H'*s;
        %t = 1./(1+lamda/(1-lamda).*((mu+u_r)./u_r).^(L/2).*exp(-diag(r*r')./(2*u_r)+diag(r*r')./(2*u_r+2*mu)));
        t = 1./(1+lamda/(1-lamda).*((mu+u_r)./u_r).^(1/2).*exp(-diag(r*r')./(2*L*u_r)+diag(r*r')./(2*L*u_r+2*L*mu)));
        x = diag(mu./(u_r+mu).*t)*r;
        u_x = (mu./(u_r+mu).*t+(1-t).*mu./(u_r.*(mu+u_r)).*diag(r*x')/L).*u_r;
%         for j = 1:M
%             p0(j,:)=H(j,:)*x0-s0(j,:)*reshape(u_p0(j,:,:),[L0,L0]);
%             p(j,:)=H(j,:)*x-s(j,:)*reshape(u_p(j,:,:),[L,L]);
%         end
%         p = H*x-u_p.*s;
    end
end
        