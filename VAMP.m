function [r1,u_r1,x1,alpha1] = VAMP(H,y,iterNum,mu,lamda,mu0)
%% H是信道矩阵，y是接收信号，iterNum是迭代次数，alpha为信号高斯分布均值，
%% mu为信号高斯分布方差，lamda为伯努利分布参数，mu0为噪声信号方差
    [M,Nc] = size(H);
    L = size(y,2);
    [U,S,V] = svd(H);
    S1 = [S;zeros(Nc-M,Nc)];
%     y_hat = U'*y;
    r1 = zeros(Nc,L) + 0.023;  %r1为初始化的x1似然变量均值
    u_r1 = ones(Nc,1);       %u_r1为r1的方差，在VAMP中所有消息的似然方差都一致，这大概是为了避免出现求x2先验时的方差为负数的情况。
    rho = 0.97;           %阻尼系数 
    x1 = zeros(Nc,L);            %x1初始值
    fi1 = inv(H'*H+1e-3*eye(Nc));
    sigma = lamda^2*1e-6+(1-lamda)^2*1;   %假设信号先验服从两个方差不同的高斯分布的线性组合，这是组合后的总体方差
    for i = 1:iterNum
        %t = 1./(1+lamda/(1-lamda)*((mu+u_r1)./u_r1).^(L/2).*exp(-mu*diag(r1*r1')./(2*u_r1.*(mu+u_r1))));
        %x1 = rho*diag(t.*mu./(mu+u_r1))*r1+(1-rho)*x1;
        %alpha1 = (mu./(mu+u_r1)).*t+(1-t).*mu./(u_r1.*(mu+u_r1)).*diag(r1*x1')/L;
        x1 = sigma*r1./(sigma+u_r1);
        alpha1 = sigma./(sigma+u_r1);
%         alpha1 = (mu./(mu+u_r1))./(1+t)+diag(diag((mu.*t./(mu+u_r1))./(abs(1+t).^2).*mu./(u_r1.*(mu+u_r1)))*(r1*r1'))/L;
%         alpha1 = diag((mu./(mu+u_r1))./(1+t))+diag((mu.*t./(mu+u_r1))./(abs(1+t).^2).*mu./(u_r1.*(mu+u_r1)))*(r1*r1')/L;
        eta1 = alpha1.*u_r1;
        u_r2 = real(eta1.*u_r1./(u_r1-eta1));
%         if(u_r2<1e-11)
%             u_r2 = 1e-11;
%         elseif(u_r2>1e11)
%             u_r2 = 1e11;
%         end
        %u_r2(find(real(u_r2)<=0))=1e-8;
        r2 = (u_r1.*x1-eta1.*r1)./(u_r1-eta1);
        %x2 = V*diag(1./((1/mu0)*diag(S1'*S1)+1./u_r2))*((1/mu0)*S'*U'*y+V'*diag(1./u_r2)*r2);
        %计算x2的过程存在问题
%         x2 = inv(1/mu0*(H'*H)+diag(1./u_r2))*(1/mu*H'*y+diag(u_r2)*r2);
        %x2 = inv(mu0*fi+diag(u_r2))*(diag(u_r2)*fi*H'*y+fi*mu0*r2);
        x2 = inv(mu0*fi1+diag(u_r2))*(diag(u_r2)*fi1*H'*y+fi1*mu0*r2);
        %V*diag(1./((1/mu0)*diag(S1'*S1)+1./u_r2))*((1/mu0)*S'*U'*y+gamma2*V'*r2);  %论文中D对角阵方程错误
        alpha2 = real(diag(inv(mu0*fi1+diag(u_r2))*mu0*fi1));
        %alpha2 = diag(inv(1/mu*(H'*H)+diag((1./u_r2))))./u_r2;
        %1/Nc*sum(gamma2./((1/mu0)*diag(S1'*S1)+gamma2));
        eta2 = alpha2.*u_r2;
        u_r1 = rho*(eta2.*u_r2)./(u_r2-eta2)+(1-rho)*u_r1;
%         if(u_r1<1e-11)
%             u_r1 = 1e-11;
%         elseif(u_r1>1e11)
%             u_r1 = 1e11;
%         end
        %u_r1(find(real(u_r1)<=0)) = 1e-8;
        r1 = (u_r2.*x2-eta2.*r2)./(u_r2-eta2);
%         r2 = (x1-alpha1*r1)./(1-alpha1);
%         u_r2 = gamma1*alpha1/(1-alpha1);
%         alpha2 = 1/Nc*sum(1./(1+mu0*u_r2*diag(S'*S)));
%         gamma1 = u_r2*alpha2/(1-alpha2);
%         r1 = r2+1/(1-alpha2)*V*inv(S'*S+mu0/u_r2*eye(Nc))*S*(y_hat-S*V'*r2);
    end
end