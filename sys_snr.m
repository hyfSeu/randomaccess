function [snr] = sys_snr(p_max, d0)
    c0 = 34.5+20*log10(d0);
    noise_figure = 9;
    n0 = -174;
    bw = 1e7;
    snr = p_max-c0-noise_figure-10*log10(bw)-n0;
    snr = 10^(snr/10);
end