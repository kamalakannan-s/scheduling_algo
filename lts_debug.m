clear all
close all
clc
% LTS
lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
lts_t = ifft([lts_f], 64);
act_coeff = complex(randn(1),randn(1));
rx = act_coeff*lts_t;
rx_f = fft(rx,64);
corr = rx_f.*[lts_f];
coeff = mean(corr);
