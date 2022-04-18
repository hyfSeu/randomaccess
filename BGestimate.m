function [Hhat,Hvar,valInH] = BGestimate(Rhat,Rvar,delta,Hmean,Hmu)
    Kappa = 1./(1+(1-delta)/delta*sqrt((Hmu+Rvar)./Rvar).*exp(-(abs(Rhat).^2)./(2*Rvar)+(abs(Rhat-Hmean).^2)./(2*(Hmu+Rvar))));%exp(-Hmu.*(abs(Rhat).^2)./(2*Rvar.*(Rvar+Hmu))))
    Hhat = Kappa.*(Hmu.*Rhat+Rvar.*Hmean)./(Hmu+Rvar);
    Hvar = ((1-Kappa)./Kappa).*(abs(Hhat).^2)+Kappa.*(Hmu.*Rvar./(Hmu+Rvar));
%     if(any(isnan(Hvar(:))))
%         stop = 1;
%     end
    temp = 0.5*(log(Rvar./(Rvar+Hmu))+Hmu./(Rvar+Hmu)-abs(Hmu.*(Rhat-Hmean)./(Rvar+Hmu)).^2./Hmu);
    valInH = Kappa.*temp+Kappa.*log(max(1e-8,delta)./max(Kappa,1e-8))+(1-Kappa).*log(max(1e-8,(1-delta))./max((1-Kappa),1e-8));
end