clear all;
close all;
clc;
M = 2;
cell_number = 6;
snr = 3;
Total_cell = 7;%neighbor
cell_ind = mod((cell_number-1),Total_cell)+1;
pow_dif = 3;
H_sv = 1/sqrt(2)*complex(randn(M,M),randn(M,M));
% H_sv = eye(M);
H_vv = (10 ^(-pow_dif/20))/sqrt(2)*complex(randn(M,M),randn(M,M));
% H_vv = eye(M);
% % % % H_vv = (10 ^(-pow_dif/20))*H_sv;
rank = 2;
[U,S,V] = svd(H_vv);
precoder_ideal = V(:,1:rank);
% H_vv = eye(2);
sym_size = 50;
MOD_ORDER = 1;%victims mod order
MOD_ORDER_int = 1; %interferer's mod order
DO_APPLY_CFO_CORRECTION = 0;
[tx_data tx_vec_air spacing tx_syms] = tx_data_gen_SIC(M,cell_ind,Total_cell,H_vv,rank,sym_size,MOD_ORDER);
%create interference
cell_number_int = 3;
cell_ind_int = mod((cell_number_int-1),Total_cell)+1;
pos = zeros(1,Total_cell);
loc = [cell_ind cell_ind_int];
pos(loc) = 1;
H_ss = 1/sqrt(2)*complex(randn(M,M),randn(M,M));
rank_int = 1;
[tx_data_int tx_vec_air_int spacing tx_syms_int] = tx_data_gen_SIC(M,cell_ind_int,Total_cell,H_ss,rank_int,sym_size,MOD_ORDER_int);

%rx vec
rx_vec_air_ref = H_vv*tx_vec_air;
noise_power = var(rx_vec_air_ref(1,:)) * 10 ^(-snr/20);
noise_mat = noise_power*complex(randn(size(rx_vec_air_ref)), randn(size(rx_vec_air_ref)));
rx_vec_air_ref = rx_vec_air_ref + noise_mat;
rx_vec_air = rx_vec_air_ref + H_sv*tx_vec_air_int;
% rx_vec_air = H_vv*tx_vec_air;


for i = 1:M
    raw_rx_dec(i,:)=downsample_SIC(rx_vec_air(i,:));
end
[mimo_training_ind, payload_ind, rx_cfo_corr_t]=packet_detection_SIC(raw_rx_dec(1,:),spacing,Total_cell,DO_APPLY_CFO_CORRECTION);
% [h_f,H] = channel_estimation_SIC(M,cell_ind,mimo_training_ind,spacing,raw_rx_dec);
[H, N0] = channel_estimation_all_SIC(M,cell_ind,mimo_training_ind,spacing,raw_rx_dec,Total_cell,pos);
%decoding
[ber_int, rx_data_int] = SIC_decoder(H(:,:,cell_ind_int),H_ss,H(:,:,cell_ind),rank_int,raw_rx_dec,payload_ind,sym_size,MOD_ORDER_int,precoder_ideal,M,tx_data_int);
%checking for successful decoding
if(ber_int == 0)
    fprintf("\n decoding is successful");
end
% [tx_data_int tx_vec_air_int spacing tx_syms_int] = tx_data_gen_SIC(M,cell_ind_int,Total_cell,rank_int,sym_size,MOD_ORDER_int,);

[syms_precoded] = datatosymprecoded_SIC(rx_data_int, rank_int, MOD_ORDER_int, H_ss);
[interferer_vec_air] = interferer_data_gen_SIC(M,cell_ind_int,Total_cell,rank_int,sym_size,MOD_ORDER_int,syms_precoded);
%cancelling interferer
rx_vec_air = rx_vec_air - H_sv*interferer_vec_air;
for i = 1:M
    raw_rx_dec(i,:)=downsample_SIC(rx_vec_air(i,:));
end
% [mimo_training_ind, payload_ind, rx_cfo_corr_t]=packet_detection_SIC(raw_rx_dec(1,:),spacing,Total_cell,DO_APPLY_CFO_CORRECTION);
% [h_f,H] = channel_estimation_SIC(M,cell_ind,mimo_training_ind,spacing,raw_rx_dec);
% [H] = channel_estimation_all_SIC(M,cell_ind,mimo_training_ind,spacing,raw_rx_dec,Total_cell);
%preparing matrix
% for i = 1:M
%     reconstruct(:,:,i) = payload_mat_gen_SIC(syms_precoded(i,:),i,sym_size);
% end
%decoding victim
[U,S,V] = svd(H(:,:,cell_ind));
precoder = V(:,1:rank);
decoder = U(:,1:rank);
[rx_syms_case] = ofdmtodataconv_SIC(raw_rx_dec,payload_ind,sym_size,M);
rx_data = decoder'*rx_syms_case;
trel = poly2trellis(7, [171 133]);
for i = 1:rank
figure(i+2);
scatter(real(rx_data(i,:)), imag(rx_data(i,:)),'filled');
title(' Signal Space of received victim bits after hard SIC');
xlabel('I'); ylabel('Q');
Demap_out_case(i,:) = demapper(rx_data(i,:),MOD_ORDER,1);
rx_data_final(i,:)= vitdec(Demap_out_case(i,:),trel,7,'trunc','hard');
[number,ber(i)] = biterr(tx_data(:,i).',rx_data_final(i,:));
end
[U_ss,S,V_ss] = svd(H_ss);
V = V_ss(:,1:rank_int);
G_sv = H(:,:,cell_ind_int)*V_ss;
G_vv = H(:,:,cell_ind)*precoder;
for m = 1:rank_int
    SINR_sv(m) = SINR_calc(G_sv,zeros(M,M),G_vv,m,M,N0);
end
for m = 1:rank
    pre_SIC_SINR_vv(m) = SINR_calc(G_vv,G_sv,zeros(M,M),m,M,N0);
    post_SIC_SINR_vv(m) = SINR_calc(G_vv,zeros(M,M),zeros(M,M),m,M,N0);
end
rate_sv = 0.8*log2(1+SINR_sv/1.333);
rate_vv = 0.8*log2(1+post_SIC_SINR_vv/1.333);