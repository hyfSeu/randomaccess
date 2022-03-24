function [Hhat,Hvar,valInH] = BGestimate(Rhat,Rvar,delta,Hmu)
    Kappa = 1./(1+(1-delta)/delta*sqrt((Hmu+Rvar)./Rvar).*exp(-Hmu.*(abs(Rhat).^2)./(2*Rvar.*(Rvar+Hmu))));
    Hhat = Kappa.*(Hmu.*Rhat)./(Hmu+Rvar);
    Hvar = Kappa.*(abs((Hmu.*Rhat)./(Hmu+Rvar)).^2+Hmu.*Rvar./(Hmu+Rvar))-abs(Hhat).^2;
    temp = 0.5*(log(Rvar./(Rvar+Hmu))+Hmu./(Rvar+Hmu)-abs(Hmu.*Rhat./(Rvar+Hmu)).^2./Hmu);
    valInH = Kappa.*temp+Kappa.*log(max(1e-8,delta)./max(Kappa,1e-8))+(1-Kappa).*log(max(1e-8,(1-delta))./max((1-Kappa),1e-8));
end