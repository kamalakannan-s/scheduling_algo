function [tx_data tx_vec_air spacing tx_syms] = tx_data_gen_SIC(M,cell_ind,Total_cell,H,rank,sym_size,MOD_ORDER)
[U,S,V] = svd(H);
precoder = V(:,1:rank);
decoder = U(:,1:rank);
pilots(:,1)= [1 1 -1 1].';
pilots(:,2)= [0 0 0 0].';
SC_IND_PILOTS           = [8 22 44 58];                           % Pilot subcarrier indices
SC_IND_DATA             = [2:7 9:21 23:27 39:43 45:57 59:64];     % Data subcarrier indices
N_SC                    = 64;                                     % Number of subcarriers
CP_LEN                  = 16; 
for i = 1:rank
    [tx_data(:,i),tx_code(:,i),preamble_legacy_A,preamble_legacy_B,preamble_mimo_A,preamble_mimo_B] = Tx_data_generation(M,rank,sym_size,MOD_ORDER);
end
for i = 1:rank
tx_sym(:,i) = mapping(tx_code(:,i)', MOD_ORDER, 1);
end
tx_syms = precoder*tx_sym.';
for i = 1:M
    [tx_vec_payload(:,i)]=generate_ofdm_payload(tx_syms(i,:),SC_IND_DATA,SC_IND_PILOTS,sym_size,pilots(:,i),N_SC,CP_LEN);
end
%interpolation filter
interp_filt2 = zeros(1,43);
interp_filt2([1 3 5 7 9 11 13 15 17 19 21]) = [12 -32 72 -140 252 -422 682 -1086 1778 -3284 10364];
interp_filt2([23 25 27 29 31 33 35 37 39 41 43]) = interp_filt2(fliplr([1 3 5 7 9 11 13 15 17 19 21]));
interp_filt2(22) = 16384;
interp_filt2 = interp_filt2./max(abs(interp_filt2));
INTERP_RATE = 2;
TX_SCALE = 1;
CFO_FLAG = 0;
%end
[preamble,spacing] = preamble_generator(M,cell_ind,Total_cell);
for i = 1:M
[tx_vec_air(:,i)]=generate_mimo(preamble(i,:),tx_vec_payload(:,i).', INTERP_RATE, TX_SCALE,CFO_FLAG,interp_filt2);
end
tx_vec_air = tx_vec_air.';
end