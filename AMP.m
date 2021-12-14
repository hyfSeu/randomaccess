% function [r,u_r,x,u_x] = AMP(H,y,iterNum,mu,lamda,mu0)
%     [M,Nc] = size(H);
%     x = zeros(Nc,1) + 0.001;
%     u_x = zeros(Nc,1)+0.09;
%     s = zeros(M,1);
%     for i =1:iterNum
%         u_p = (abs(H).^2) * u_x;
%         p = H * x - u_p.*s;
%         z = H * x;
%         s = (y - p)./(mu0 + u_p);
%         u_s = 1./(mu0 + u_p);
%         u_r = 1./((abs(H).^2)' * u_s);
%         r = x + u_r.*(H'*s);
%         t = 1./(1+((lamda/(1-lamda)).*sqrt((u_r+mu)./u_r).*exp(-(abs(r).^2)./(2*u_r)+(abs(r).^2)./(2*(u_r+mu)))));
%         x = ((mu*r)./(mu+u_r)).*t; %((mu*r)./(mu+u_r))./(1+((lamda/(1-lamda)).*sqrt((u_r+mu)./u_r).*exp(-(abs(r).^2)./(2*u_r)+(abs(r).^2)./(2*(u_r+mu)))));
%         u_x = ((mu*u_r)./(mu+u_r)).*t+ ((mu*u_r)./(mu+u_r)).*t.*(1-t).*((abs(r).^2)./u_r-(abs(r).^2)./(u_r+mu)); %mu*u_r./(alpha*u_r+r*mu).*x+((mu*r+u_r*alpha)./(mu+u_r)-x).*x;
%     end
% end


function [r,u_r,x,sigma] = AMP(H,y,iterNum,mu,lamda,mu0)
%% H是信道矩阵，y是接收信号，iterNum是迭代次数，alpha为信号高斯分布均值，
%% mu为信号高斯分布方差，lamda为伯努利分布参数，mu0为噪声信号方差
    [M,Nc] = size(H);
    L = size(y,2);
    if(L<2)
        L0 = L;
    else
        L0 = 2;
    end
    K = L/L0;
    x = zeros(Nc,L) + 0.023;
    x0 = x(:,1:L0);
    u_x = zeros(Nc,L,L);
    u_x0 = u_x(:,1:L0,1:L0);
    for i =1:Nc
        u_x0(i,:,:) = 0.09*eye(L0);
    end
    s0 = zeros(M,L0);
    s = zeros(M,L);
    u_s0 = zeros(M,L0,L0);
    u_s = zeros(M,L,L);
    r0 = zeros(Nc,L0);
    r = zeros(Nc,L);
    u_r0 = zeros(Nc,L0,L0);
    u_r = zeros(Nc,L,L);
    p0 = zeros(M,L0);
    p = zeros(M,L);
    sigma = zeros(L,Nc,Nc);
    for i = 1:iterNum
        u_p0 = pagemtimes((abs(H).^2),u_x0);
        u_p = pagemtimes((abs(H).^2),u_x);
        for j = 1:M
            p0(j,:)=H(j,:)*x0-s0(j,:)*reshape(u_p0(j,:,:),[L0,L0]);
            p(j,:)=H(j,:)*x-s(j,:)*reshape(u_p(j,:,:),[L,L]);
        end
%         p = H*x-u_p.*s;
        z = H*x0;
        for j = 1:M
            u_s0(j,:,:)=inv(mu0*eye(L0)+reshape(u_p0(j,:,:),[L0,L0]));
            u_s(j,:,:)=inv(mu0*eye(L)+reshape(u_p(j,:,:),[L,L]));
            s0(j,:)=(y(j,1:L0)-p0(j,:))*reshape(u_s0(j,:,:),[L0,L0]);
            s(j,:)=(y(j,1:L)-p(j,:))*reshape(u_s(j,:,:),[L,L]);
        end
%         s = (y-p)./(u_p+mu0);
%         u_s = 1./(u_p+mu0);
        u_r0 = pagemtimes((abs(H).^2)',u_s0);
        u_r = pagemtimes((abs(H).^2)',u_s);
        for j = 1:Nc
            u_r0(j,:,:) = inv(reshape(u_r0(j,:,:),[L0,L0]));
            u_r(j,:,:) = inv(reshape(u_r(j,:,:),[L,L]));
            r0(j,:) = x0(j,:)+H(:,j)'*s0*reshape(u_r0(j,:,:),[L0,L0]);
            r(j,:) = x(j,:)+H(:,j)'*s*reshape(u_r(j,:,:),[L,L]);
        end
%         u_r = 1./((H.*((H.')'))'*u_s);
%         r = x+u_r.*(H'*s);
        for j = 1:Nc  %%这里修改了u_x的定义，r_H*r放在前面
            %在L较大的情况下，会出现行列式无穷大的情况，需要增加一个判别分别处理；修改点：将flag从det(inv(reshape(u_r(j,:,:),[L,L])))
            %改为det((mu*eye(L)+reshape(u_r(j,:,:),[L,L]))*inv(reshape(u_r(j,:,:),[L,L])))
            flag = det((mu*eye(L0)+reshape(u_r0(j,:,:),[L0,L0]))*inv(reshape(u_r0(j,:,:),[L0,L0])));
            if(abs(flag)==inf)
                t = 1e-25;
%             elseif flag==inf
%                 t = mu/(1+lamda/(1-lamda)*sqrt(det(mu*eye(L)))*exp(-0.5*r(j,:)*inv(reshape(u_r(j,:,:),[L,L]))*r(j,:)'+0.5*r(j,:)*inv(mu*eye(L)+reshape(u_r(j,:,:),[L,L]))*r(j,:)'));
            else
                t = mu/(1+lamda/(1-lamda)*sqrt(flag)*exp(-0.5*r0(j,:)*inv(reshape(u_r0(j,:,:),[L0,L0]))*r0(j,:)'+0.5*r0(j,:)*inv(mu*eye(L0)+reshape(u_r0(j,:,:),[L0,L0]))*r0(j,:)'));
            end
            x0(j,:) = r0(j,:)*t*inv(mu*eye(L0)+reshape(u_r0(j,:,:),[L0,L0]));
            x(j,:) = r(j,:)*t*inv(mu*eye(L)+reshape(u_r(j,:,:),[L,L]));
            u_x0(j,:,:) = (t*inv(mu*eye(L0)+reshape(u_r0(j,:,:),[L0,L0]))+x0(j,:)'*r0(j,:)*(1-t/mu)*(inv(reshape(u_r0(j,:,:),[L0,L0]))-inv(mu*eye(L0)+reshape(u_r0(j,:,:),[L0,L0]))))*reshape(u_r0(j,:,:),[L0,L0]);
            u_x(j,:,:) = kron(eye(K),reshape(u_x0(j,:,:),[L0,L0]));
        end
    end
    for i = 1:L
        sigma(i,:,:) = diag(reshape(u_x(:,i,i),[Nc,1]));
    end
end
        