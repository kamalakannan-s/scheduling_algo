function [success,N0] = SIC_PHY(rank,rank_int,mod_order,mod_order_int,cell_number,cell_number_int,H_ss,H_vv,H_sv)
M = 2;
MOD_ORDER = mod_order;
MOD_ORDER_int = mod_order_int;
snr = 10;%to be removed
Total_cell = 7;%neighbor
cell_ind = mod((cell_number-1),Total_cell)+1;
[U,S,V] = svd(H_vv);
precoder_ideal = V(:,1:rank);
% H_vv = eye(2);
sym_size = 50;
DO_APPLY_CFO_CORRECTION = 0;
[tx_data tx_vec_air spacing tx_syms] = tx_data_gen_SIC(M,cell_ind,Total_cell,H_vv,rank,sym_size,MOD_ORDER);
%create interference
cell_ind_int = mod((cell_number_int-1),Total_cell)+1;
pos = zeros(1,Total_cell);
loc = [cell_ind cell_ind_int];
pos(loc) = 1;
[tx_data_int tx_vec_air_int spacing tx_syms_int] = tx_data_gen_SIC(M,cell_ind_int,Total_cell,H_ss,rank_int,sym_size,MOD_ORDER_int);

%rx vec
rx_vec_air_ref = H_vv*tx_vec_air;
noise_power = var(rx_vec_air_ref(1,:)) * 10 ^(-snr/20);
noise_mat = noise_power*complex(randn(size(rx_vec_air_ref)), randn(size(rx_vec_air_ref)));
rx_vec_air_ref = rx_vec_air_ref;% + noise_mat;
rx_vec_air = rx_vec_air_ref + H_sv*tx_vec_air_int;
% rx_vec_air = H_vv*tx_vec_air;


for i = 1:M
    raw_rx_dec(i,:)=downsample_SIC(rx_vec_air(i,:));
end
[mimo_training_ind, payload_ind, rx_cfo_corr_t, fail]=packet_detection_SIC(raw_rx_dec(1,:),spacing,Total_cell,DO_APPLY_CFO_CORRECTION);
% [h_f,H] = channel_estimation_SIC(M,cell_ind,mimo_training_ind,spacing,raw_rx_dec);
[H, N0] = channel_estimation_all_SIC(M,cell_ind,mimo_training_ind,spacing,raw_rx_dec,Total_cell,pos);
%decoding
[ber_int, rx_data_int] = SIC_decoder(H(:,:,cell_ind_int),H_ss,H(:,:,cell_ind),rank_int,raw_rx_dec,payload_ind,sym_size,MOD_ORDER_int,precoder_ideal,M,tx_data_int);
%checking for successful decoding
if(ber_int == 0)
    fprintf("\n decoding is successful");
end


[syms_precoded] = datatosymprecoded_SIC(rx_data_int, rank_int, MOD_ORDER_int, H_ss);
[interferer_vec_air] = interferer_data_gen_SIC(M,cell_ind_int,Total_cell,rank_int,sym_size,MOD_ORDER_int,syms_precoded);
%cancelling interferer
rx_vec_air_1 = rx_vec_air - H_sv*interferer_vec_air;
for i = 1:M
    raw_rx_dec_1(i,:)=downsample_SIC(rx_vec_air_1(i,:));
end
%decoding victim
[U,S,V] = svd(H(:,:,cell_ind));
precoder = V(:,1:rank);
decoder = U(:,1:rank);
[rx_syms_case] = ofdmtodataconv_SIC(raw_rx_dec_1,payload_ind,sym_size,M);
rx_data = decoder'*rx_syms_case;
trel = poly2trellis(7, [171 133]);
for i = 1:rank
figure(i+2);
scatter(real(rx_data(i,:)), imag(rx_data(i,:)),'filled');
title(' Signal Space of received victim bits after hard SIC');
xlabel('I'); ylabel('Q');
Demap_out_case(i,:) = demapper(rx_data(i,:),MOD_ORDER,1);
rx_data_final(i,:)= vitdec(Demap_out_case(i,:),trel,7,'trunc','hard');
[~,ber(i)] = biterr(tx_data(:,i).',rx_data_final(i,:));
end
BER = [ber_int ber];
if sum(BER)==0
    success = 1;
else
    success = 0;
end
N0 = 0.2;
end