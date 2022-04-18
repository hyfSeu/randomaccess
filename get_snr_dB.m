function snr = get_snr_dB(d,d0,pmax,pl_exp)
%% this function is used to compute the large-scale fading factor with normalized
%% noise variance, as a snr dB form
    if(~isempty(find(d>300,1)) || ~isempty(find(d<0,1)))
        error("the distance is wrong!");
    end
    c0 = 34.5+20*log10(d0);
    noise_figure = 9;
    n0 = -174;
    bw = 1e7;
    snr0 = pmax-c0-noise_figure-10*log10(bw)-n0;
    snr = 10*log10(2./(1+(1+d/d0).^pl_exp)) + snr0;
end