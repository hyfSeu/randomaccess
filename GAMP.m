clear;
Nc = 512;
Ka = 13;
M = 120;
sigma = 1;
mu = 1;
nu = 1;
one_indx = randperm(Nc,Ka);
one_indx = sort(one_indx);
message = zeros(Nc,1);
message(one_indx) = normrnd(mu,nu,Ka,1);
signals = 1-2*message;
%y = awgn(signals,4,'measured');
%A = PartHadamardMtx(M,Nc);
A = sqrt(100)*normrnd(0,1,M,Nc);
y = A * message + normrnd(0,sigma,M,1);
lamda = (Nc - Ka)/Nc;

x = zeros(Nc,1)+0.001;
u_x = zeros(Nc,1)+0.009;
s = zeros(M,1);
iterNums = 500;
for i =1:iterNums
    u_p = (A .* A) * u_x;
    p = A * x - u_p.*s;
    z = A * x;
    s = (y - p)./(sigma^2 + u_p);
    u_s = 1./(sigma^2 + u_p);
    u_r = 1./((A.* A)' * u_s);
    r = x + u_r.*(A'*s);
    %x = 1./((lamda/(1-lamda)).*exp((1-2*r)./(2*u_r))+1);
    %u_x = x.*(1-x);
    x = ((nu^2*r+u_r*mu)./(nu^2+u_r))./(1+((lamda/(1-lamda)).*sqrt((u_r+nu^2)./u_r).*exp(-(r.*r)./(2*u_r)+(mu-r).*(mu-r)./(2*(u_r+nu^2)))));
    u_x = nu^2*u_r./(mu*u_r+r*nu^2).*x+((nu^2*r+u_r*mu)./(nu^2+u_r)-x).*x;
    %temp1 = (r.*r)./(2 * u_r) - (r - 1).*(r - 1)./(2*(0.01+u_r))-9.58;
    %temp1 = (2*r-1)./(2*u_r)-4.06;
    %index = find(temp1<0);
    %x(index) = 1;
    %u_x(index) = 0;
end
t = message(one_indx);
tt = x(one_indx);