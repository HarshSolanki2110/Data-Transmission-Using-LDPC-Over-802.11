function msg_cap = BPSK_Calcu(EbN0,msg)
    R = 1;
    N = size(msg,1);
    sigma = sqrt(1/ (2*R*EbN0));
    s = 1 - 2 * msg;
    msg_cap = s + sigma .* randn(N,1);
end