function[h_f,H]=channel_estimation_SIC(M,cell_ind,mimo_training_ind,spacing,rx_data)

% LTS
lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];


SC_IND_DATA             = [2:7 9:21 23:27 39:43 45:57 59:64]; 
for i = 1:M
    lts_ind_start(i) = mimo_training_ind + 32 + spacing*(cell_ind-1) + 96*(i-1);
    lts_ind_end(i) = lts_ind_start(i) + 64 - 1;
end
for i = 1:M^2
    h((floor((i-1)/M))+1,mod((i-1),M)+1,:) = rx_data((floor((i-1)/M))+1,lts_ind_start(mod((i-1),M)+1):lts_ind_end(mod((i-1),M)+1));
    ref((floor((i-1)/M))+1,mod((i-1),M)+1,:) = lts_f;
end
h_f = fft(h,64,3);
h_corr = h_f.*ref;
h_est = h_corr(:,:,SC_IND_DATA);
H = mean(h_est,3);
end