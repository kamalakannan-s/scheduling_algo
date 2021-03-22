function[H, N0]=channel_estimation_all_SIC(M,cell_ind,mimo_training_ind,spacing,rx_data,Total_cell,pos)

%find location for noise measurement
free_slot = find(~pos);
occupied = length(free_slot);
if (occupied == 1)
    link_class = 0;
else
    link_class = 1;
end
N_slot = free_slot(1);
% LTS
lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];


SC_IND_DATA             = [2:7 9:21 23:27 39:43 45:57 59:64]; 
for j = 1:Total_cell
for i = 1:M
    lts_ind_start(i) = mimo_training_ind + 32 + spacing*(j-1) + 96*(i-1);
    lts_ind_end(i) = lts_ind_start(i) + 64 - 1;
end
for i = 1:M^2
    h((floor((i-1)/M))+1,mod((i-1),M)+1,:) = rx_data((floor((i-1)/M))+1,lts_ind_start(mod((i-1),M)+1):lts_ind_end(mod((i-1),M)+1));
    ref((floor((i-1)/M))+1,mod((i-1),M)+1,:) = lts_f;
end
%noise power measurement
noise_ind_start = mimo_training_ind + 32 + spacing*(N_slot-1);
noise_ind_end = noise_ind_start + 64 - 1;
N0 = var(rx_data(1,noise_ind_start:noise_ind_end));
%end
h_f = fft(h,64,3);
h_corr = h_f.*ref;
h_est = h_corr(:,:,SC_IND_DATA);
H(:,:,j) = mean(h_est,3);
end
end