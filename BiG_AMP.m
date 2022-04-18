function [resX,resXvar,resH,resHvar] = BiG_AMP(Y,noise_power,K,delta,Xmean,Xmu,Hmean,Hmu,init_X,init_Xvar,init_H,init_Hvar)
%% this function is used to perform BiG-AMP. We assume the X obey gaussian distribution with 0 and Xmu
%% and H obey BG distribution with factor delta,0,Hmu
    resX = [];
    resXvar = [];
    resH = [];
    resHvar = [];
    [L,M] = size(Y);
    tao = 1e-8;                                      %迭代终止阈值
    IterNum = 8000;                                  %AMP算法迭代次数
    IterOut = 30;
    IterMin = 30;
    stepWindow = 1;                                 %size of stepWindw
    valOpt = [];                                     %set of succeed step results
    stepTol = -1;
    stepMin = 0.05;
    stepMax = 0.45;
    adaptStep = true;
    stepIncr = 1.1;
    stepDecr = 0.5;
    maxStepDecr = 0.5;
    PvarMin = 1e-13;                                    %% warning!!!
    ZvarToPvarMax = 0.99;                              %% warning!!!
    varThresh = 1e6;                                   %% warning!!!
    HvarMin = 0;                                    %% warning!!!
    XvarMin = 0;                                    %% warning!!!
    cfactorMin = inf;
    %% 循环迭代
    % initialize Xhat,Hhat,Xvar,Hvar,Shat
    %Lamda = repmat(lamda,1,M);                        %var of channels
    for m = 1:IterOut
        Pvar_mean = 1;
        valIn = 0;
        val = zeros(IterNum,1);
        failCount = 0;
        step = 1;
        beta = 0.05;                                        %damping系数

        %X variables
        if isempty(init_X)
            Xhat = Xmean+sqrt(Xmu).*(randn(L,K)+1j*randn(L,K))/sqrt(2);
        else
            Xhat = init_X;
        end
        if isempty(init_Xvar)
            Xvar = 10*Xmu;
        else
            Xvar = init_Xvar;
        end
        XhatBar = zeros(L,K);
        XhatOpt = zeros(L,K);
        XhatBarOpt = zeros(L,K);

        %H variables
        if isempty(init_H)
            Hhat = delta*Hmean;
        else
            Hhat = init_H;
        end
        if isempty(init_Hvar)
            Hvar = 10*Hmu;
        else
            Hvar = init_Hvar;
        end
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

           if(it>1 && adaptStep)
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
           [Hhat,Hvar,valInH] = BGestimate(Rhat,Rvar,delta,Hmean,Hmu);

           %caclulate the x variables
           [Xhat,Xvar,valInX] = Gestimate(Qhat,Qvar,Xmean,Xmu);

           %update valIn
           valIn = sum(valInH(:))+sum(valInX(:));

           %don't stop before minimum iteration count
           if it<IterMin
               stop = false;
           end
        end
    %     correct_col = sqrt(Xmu(1,:)).*(ones(1,K)+1j*ones(1,K))/sqrt(2);
    %     Q = find_permutation(correct_col,Xhat(1,:));
    %     Qinv = 1./Q;
    %     Xhat = Xhat*diag(Qinv);
    %     Hhat = diag(Q)*Hhat;
        Hhat(abs(Hhat)<1e-2)=0;
        Hhats = abs(Hhat);
        Hhats(Hhats<0.25) = 0;
        Hhats(Hhats>0.75)=1;
        check_mtx = [];
        uni_code = [];
        if(delta == 0.25)
            P = M/4;
            for i = 1:P
                temp_code = zeros(K,4);
                check_mtx = [check_mtx,sum(Hhats(:,4*i-3:4*i),2)];
                %[biggest_cidx,biggest_ridx]=find(abs(Hhat(:,4*i-3:4*i)) == max(abs(Hhat(:,4*i-3:4*i)),[],2));
                temp_code(abs(Hhat(:,4*i-3:4*i)) == max(abs(Hhat(:,4*i-3:4*i)),[],2))=1;
                uni_code = [uni_code,temp_code];
            end
        elseif delta == 0.125
            P = M/8;
            for i = 1:P
                temp_code = zeros(K,8);
                check_mtx = [check_mtx,sum(Hhats(:,8*i-7:8*i),2)];
                temp_code(abs(Hhat(:,8*i-7:8*i)) == max(abs(Hhat(:,8*i-7:8*i)),[],2))=1;
                uni_code = [uni_code,temp_code];
            end
        else
            error("wrong delta!");
        end
        recheck = pdist2(ones(1,P),check_mtx,'hamming');
        cfactor = sum(recheck);
        [UEnum,~] = size(unique(uni_code,'rows'));
        if (cfactor<cfactorMin) && (UEnum == K)
            cfactorMin = cfactor;
            resH = Hhat;
            resX = Xhat;
            resHvar = Hvar;
            resXvar = Xvar;
        elseif (isempty(resH)) && (m==IterOut)
            resH = Hhat;
            resX = Xhat;
            resHvar = Hvar;
            resXvar = Xvar;
        end
        if (cfactorMin<0.5) && (UEnum == K)
            break;
        end
    end
end