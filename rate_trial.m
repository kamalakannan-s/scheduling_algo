clear all;
close all;
clc;
pow_dif = 5;
M = 2;
% H1 = 1/sqrt(2)*complex(randn(M,M),randn(M,M));
H1 = eye(M);
H2 = (10 ^(-pow_dif/20))/sqrt(2)*complex(randn(M,M),randn(M,M));
H2 = zeros(M,M);
H3 = (10 ^(-(pow_dif+3)/20))/sqrt(2)*complex(randn(M,M),randn(M,M));
H3 = zeros(M,M);

sinr(1) = SINR_calc(H1,H2,H3,1,M)
sinr_db(1) = 10*log10(sinr(1))
sinr(2) = SINR_calc(H2,H1,H3,1,M)
sinr_db(2) = 10*log10(sinr(2))
sinr(3) = SINR_calc(H3,H1,H2,1,M)
sinr_db(3) = 10*log10(sinr(3))
rate = 0.8*log2(1+sinr/1.333)