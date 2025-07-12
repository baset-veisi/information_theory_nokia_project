function c_sh = ShannonLimitAttenuated(SNR_vec,alpha)
c_sh = zeros(length(SNR_vec),1);
for i = 1:length(SNR_vec)
    snr = 10^(SNR_vec(i)/10);
    if (1-alpha^2)/alpha^2 > snr
        c_sh(i) = log2(1+snr);
    else
        beta = (1/alpha^2 - 1)/(2*snr) + 0.5;
        c_sh(i) = log2((1+beta*snr)*(1+alpha^2*(1-beta)*snr));
    end
end
