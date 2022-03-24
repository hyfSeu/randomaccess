function [Xhat,Xvar,valInX] = Gestimate(Qhat,Qvar,Xmu)
    Xhat = (Qhat.*Xmu)./(Xmu+Qvar);
    Xvar = (Qvar.*Xmu)./(Xmu+Qvar);
    temp = Qvar./(Qvar+Xmu);
    valInX = log(temp)+(1-temp)-abs(Xhat).^2./Xmu;
end