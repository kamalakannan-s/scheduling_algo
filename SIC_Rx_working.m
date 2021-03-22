clear all;
close all;
clc;
M = 2;
iter = 1;%monte carlo iteration
snr_iter = 1;%snr steps used for initializing BER matrix
snr = 20;%snr starting dB 

pow_dif = 10;
ber = zeros(iter,snr_iter,4);
ber1 = zeros(iter,snr_iter,4);
mul_tx = 1;% 0 for simo and 1 for multiple Tx stream mimo case
modulations = ["BPSK" "QPSK" "16-QAM" "64-QAM"];
mod_i = 1;
mod_order_s = mod_i;
    if(mod_i == 4)
        offset = 2;
    else
        offset = 0;
    end
    MOD_ORDER               =  2^(mod_i-1)-offset;
    %generating interferer air data
    project_SIC_Tx_gen;
    preamble_A = [preamble_legacy_A, preamble_mimo_A zeros(1,length(preamble_mimo_A))];
    preamble_B = [preamble_legacy_B, preamble_mimo_B zeros(1,length(preamble_mimo_B))];
    tx_interferer_data_a = tx_data_a;
    tx_interferer_data_b = tx_data_b;
    [tx_vec_air_A]=generate_mimo(preamble_A,tx_payload_vec_A, INTERP_RATE, TX_SCALE,CFO_FLAG,interp_filt2);
    [tx_vec_air_B]=generate_mimo(preamble_B,tx_payload_vec_B, INTERP_RATE, TX_SCALE,CFO_FLAG,interp_filt2);
    tx_int_A = tx_vec_air_A;
    tx_int_B = tx_vec_air_B;
    %generating desire signal air data
    MOD_ORDER               =  1;
    project_SIC_Tx_gen;
    preamble_victim_A = [preamble_legacy_A, zeros(1,length(preamble_mimo_A)) preamble_mimo_A];
    preamble_victim_B = [preamble_legacy_B, zeros(1,length(preamble_mimo_B)) preamble_mimo_B];
    tx_victim_data_a = tx_data_a;
    tx_victim_data_b = tx_data_b;
%     [tx_vec_air_A]=generate_mimo(preamble_victim_A,zeros(1,length(tx_payload_vec_A)), INTERP_RATE, TX_SCALE,CFO_FLAG,interp_filt2);
%     [tx_vec_air_B]=generate_mimo(preamble_victim_B,zeros(1,length(tx_payload_vec_B)), INTERP_RATE, TX_SCALE,CFO_FLAG,interp_filt2);
    [tx_vec_air_A]=generate_mimo(preamble_victim_A,tx_payload_vec_A, INTERP_RATE, TX_SCALE,CFO_FLAG,interp_filt2);
    [tx_vec_air_B]=generate_mimo(preamble_victim_B,tx_payload_vec_B, INTERP_RATE, TX_SCALE,CFO_FLAG,interp_filt2);
    tx_desire_A = tx_vec_air_A;
    tx_desire_B = tx_vec_air_B;
    %generating channel matrix for both the interferer and signal with
    %power difference
    H1 = 1/sqrt(2)*complex(randn(M,M),randn(M,M));
    H2 = (10 ^(-pow_dif/20))/sqrt(2)*complex(randn(M,M),randn(M,M));
    %Rx at the victim
    %for sending zeros in payload
%     rx_vec_air = H1*[tx_int_A;tx_int_B]+H2*[zeros(1,length(tx_desire_A));zeros(1,length(tx_desire_B))];
    %for sending data in payload of victim
    rx_vec_air = H1*[tx_int_A;tx_int_B]+H2*[tx_desire_A;tx_desire_B];
    rx_vec_air_A = rx_vec_air(1,:);
    rx_vec_air_B = rx_vec_air(2,:);
    noise_power = var(rx_vec_air_A) * 10 ^(-snr/20);
    %mute the following 2 line for avoiding noise
    rx_vec_air_A = rx_vec_air_A +noise_power*complex(randn(1,length(rx_vec_air_A)), randn(1,length(rx_vec_air_A)));
    rx_vec_air_B = rx_vec_air_B +noise_power*complex(randn(1,length(rx_vec_air_B)), randn(1,length(rx_vec_air_B)));
    raw_rx_dec_A=generate_downsample_rx(rx_vec_air_A,INTERP_RATE, interp_filt2);
    raw_rx_dec_B=generate_downsample_rx(rx_vec_air_B,INTERP_RATE, interp_filt2);
    
%% Correlate for LTS
LTS_CORR_THRESH=.8;
DO_APPLY_CFO_CORRECTION=0;
DO_APPLY_SFO_CORRECTION=0;
DO_APPLY_PHASE_ERR_CORRECTION=0;

% For simplicity, we'll only use RFA for LTS correlation and peak
% discovery. A straightforward addition would be to repeat this process for
% RFB and combine the results for detection diversity.

% Complex cross correlation of Rx waveform with time-domain LTS
lts_corr = abs(conv(conj(fliplr(lts_t)), sign(raw_rx_dec_A)));
% % lts_corr_simo = abs(conv(conj(fliplr(lts_t)), sign(raw_rx_dec_A_SIMO)));
% Skip early and late samples - avoids occasional false positives from pre-AGC samples
lts_corr = lts_corr(32:end-32);

% Find all correlation peaks
lts_peaks = find(lts_corr > LTS_CORR_THRESH*max(lts_corr));

% Select best candidate correlation peak as LTS-payload boundary
% In this MIMO example, we actually have 3 LTS symbols sent in a row.
% The first two are sent by RFA on the TX node and the last one was sent
% by RFB. We will actually look for the separation between the first and the
% last for synchronizing our starting index.

[LTS1, LTS2] = meshgrid(lts_peaks,lts_peaks);
[lts_last_peak_index,y] = find(LTS2-LTS1 == length(lts_t));

% Stop if no valid correlation peak was found
if(isempty(lts_last_peak_index))
    fprintf('No LTS Correlation Peaks Found!\n');
%     return;
end

% Set the sample indices of the payload symbols and preamble
% The "+32" here corresponds to the 32-sample cyclic prefix on the preamble LTS
% The "+192" corresponds to the length of the extra training symbols for MIMO channel estimation
mimo_training_ind = lts_peaks(max(lts_last_peak_index)) + 32;
payload_ind = mimo_training_ind + 2*192;

% Subtract of 2 full LTS sequences and one cyclic prefixes
% The "-160" corresponds to the length of the preamble LTS (2.5 copies of 64-sample LTS)
lts_ind = mimo_training_ind-160;

if(DO_APPLY_CFO_CORRECTION)
    %Extract LTS (not yet CFO corrected)
    rx_lts = raw_rx_dec_A(lts_ind : lts_ind+159); %Extract the first two LTS for CFO
    rx_lts1 = rx_lts(-64 + [97:160]);
    rx_lts2 = rx_lts( [97:160]);

    %Calculate coarse CFO est
    rx_cfo_est_lts = mean(unwrap(angle(rx_lts2 .* conj(rx_lts1))));
    rx_cfo_est_lts = rx_cfo_est_lts/(2*pi*64);
else
    rx_cfo_est_lts = 0;
end

% Apply CFO correction to raw Rx waveforms
rx_cfo_corr_t = exp(-1i*2*pi*rx_cfo_est_lts*[0:length(raw_rx_dec_A)-1]);
rx_dec_cfo_corr_A = raw_rx_dec_A .* rx_cfo_corr_t;
rx_dec_cfo_corr_B = raw_rx_dec_B .* rx_cfo_corr_t;


% MIMO Channel Estimatation for interferer
lts_ind_TXA_start = mimo_training_ind + 32 ;
lts_ind_TXA_end = lts_ind_TXA_start + 64 - 1;

lts_ind_TXB_start = mimo_training_ind + 32 + 64 + 32 ;
lts_ind_TXB_end = lts_ind_TXB_start + 64 - 1;

rx_lts_AA = rx_dec_cfo_corr_A( lts_ind_TXA_start:lts_ind_TXA_end );
rx_lts_BA = rx_dec_cfo_corr_A( lts_ind_TXB_start:lts_ind_TXB_end );

rx_lts_AB = rx_dec_cfo_corr_B( lts_ind_TXA_start:lts_ind_TXA_end );
rx_lts_BB = rx_dec_cfo_corr_B( lts_ind_TXB_start:lts_ind_TXB_end );

rx_lts_AA_f = fft(rx_lts_AA);
rx_lts_BA_f = fft(rx_lts_BA);

rx_lts_AB_f = fft(rx_lts_AB);
rx_lts_BB_f = fft(rx_lts_BB);

%% Perform Channel estimation 
h11_sv = rx_lts_AA_f.*lts_f;
h12_sv = rx_lts_BA_f.*lts_f;
h21_sv = rx_lts_AB_f.*lts_f;
h22_sv = rx_lts_BB_f.*lts_f;
H_sv = zeros(M,M,64);
mul_sv = zeros(M,M,64);
for idx = 1:64

      H_sv(:,:,idx) = [h11_sv(idx) h12_sv(idx);h21_sv(idx) h22_sv(idx)];
      mul_sv(:,:,idx) = inv(noise_power*eye(2,2)+H_sv(:,:,idx)'*H_sv(:,:,idx))*H_sv(:,:,idx)';%MMSE combiner

end

%MIMO channel estimation for victim
lts_ind_TXA_start = mimo_training_ind + 32 + 192 ;
lts_ind_TXA_end = lts_ind_TXA_start + 64 - 1;

lts_ind_TXB_start = mimo_training_ind + 32 + 64 + 32 + 192 ;
lts_ind_TXB_end = lts_ind_TXB_start + 64 - 1;

rx_lts_AA = rx_dec_cfo_corr_A( lts_ind_TXA_start:lts_ind_TXA_end );
rx_lts_BA = rx_dec_cfo_corr_A( lts_ind_TXB_start:lts_ind_TXB_end );

rx_lts_AB = rx_dec_cfo_corr_B( lts_ind_TXA_start:lts_ind_TXA_end );
rx_lts_BB = rx_dec_cfo_corr_B( lts_ind_TXB_start:lts_ind_TXB_end );

rx_lts_AA_f = fft(rx_lts_AA);
rx_lts_BA_f = fft(rx_lts_BA);

rx_lts_AB_f = fft(rx_lts_AB);
rx_lts_BB_f = fft(rx_lts_BB);

%% Perform Channel estimation 
h11 = rx_lts_AA_f.*lts_f;
h12 = rx_lts_BA_f.*lts_f;
h21 = rx_lts_AB_f.*lts_f;
h22 = rx_lts_BB_f.*lts_f;
H_vv = zeros(M,M,64);
mul_vv = zeros(M,M,64);
for idx = 1:64

      H_vv(:,:,idx) = [h11(idx) h12(idx);h21(idx) h22(idx)];
      mul_vv(:,:,idx) = inv(noise_power*eye(2,2)+H_vv(:,:,idx)'*H_vv(:,:,idx))*H_vv(:,:,idx)';%MMSE combiner

end

%% Rx payload processing, Perform combining for 1X4 and 2X2 separately  

% Extract the payload samples (integral number of OFDM symbols following preamble)
rxA_preamble_removed = rx_dec_cfo_corr_A(payload_ind+(0:(CP_LEN+N_SC)*(N_OFDM_SYMS)-1));
rxB_preamble_removed = rx_dec_cfo_corr_B(payload_ind+(0:(CP_LEN+N_SC)*(N_OFDM_SYMS)-1));
payload_mat_A = reshape(rxA_preamble_removed,CP_LEN+N_SC,N_OFDM_SYMS);
payload_mat_B = reshape(rxB_preamble_removed,CP_LEN+N_SC,N_OFDM_SYMS);
% Remove the cyclic prefix
payload_mat_noCP_A = payload_mat_A(CP_LEN+[1:N_SC], :);
payload_mat_noCP_B = payload_mat_B(CP_LEN+[1:N_SC], :);

% Take the FFT
syms_f_mat_A = fft(payload_mat_noCP_A, N_SC, 1);
syms_f_mat_B = fft(payload_mat_noCP_B, N_SC, 1);


%*This is optional -- SFO correction*



if DO_APPLY_SFO_CORRECTION
   

else
	% Define an empty SFO correction matrix (used by plotting code below)
    pilot_phase_sfo_corr = zeros(N_SC, N_OFDM_SYMS);
end

%*This is optional* 
% Extract the pilots and calculate per-symbol phase error
if DO_APPLY_PHASE_ERR_CORRECTION
    pilots_f_mat = syms_eq_mat_pilots(SC_IND_PILOTS, :);
    pilot_phase_err = angle(mean(pilots_f_mat.*pilots_mat_A));
else
	% Define an empty phase correction vector (used by plotting code below)
    pilot_phase_err = zeros(1, N_OFDM_SYMS);
end
pilot_phase_corr = repmat(exp(-1i*pilot_phase_err), N_SC, 1);

% Apply pilot phase correction to both received streams
syms_f_mat_pc_A = syms_f_mat_A .* pilot_phase_corr;
syms_f_mat_pc_B = syms_f_mat_B .* pilot_phase_corr;

% Perform combining for MIMO 1X4 and 2X2 
% you need to apply the MIMO equalization to each subcarrier separately and then perform combining
%equalization

syms_eq_mat_A = zeros(N_SC,N_OFDM_SYMS);
syms_eq_mat_B = zeros(N_SC,N_OFDM_SYMS);
for i = 1:N_OFDM_SYMS
    for j = 1:N_SC
         correction_fac = mul_sv(:,:,j)*[syms_f_mat_pc_A(j,i) syms_f_mat_pc_B(j,i)].';
         syms_eq_mat_A(j,i) = correction_fac(1);
         syms_eq_mat_B(j,i) = correction_fac(2);
    end
end


payload_syms_mat_A = syms_eq_mat_A(SC_IND_DATA, :);
payload_syms_mat_B = syms_eq_mat_B(SC_IND_DATA, :);



%% perform demodulate or demapping post combined symbols 
rx_syms_case_1 = reshape(payload_syms_mat_A,1,numel(payload_syms_mat_A));
rx_syms_case_2 = reshape(payload_syms_mat_B,1,numel(payload_syms_mat_B));
% plot the demodulated output rx_syms_case_1 and rx_syms_case_2
figure(4);subplot(2,2,1);
scatter(real(rx_syms_case_1), imag(rx_syms_case_1),'filled');
title(' Signal Space of received strong interferer bits');
xlabel('I'); ylabel('Q');

figure(4);subplot(2,2,2);
scatter(real(rx_syms_case_2), imag(rx_syms_case_2),'filled');
title(' Signal Space of received strong interferer bits');
xlabel('I'); ylabel('Q');


% FEC decoder for the rx_syms_case_1 and rx_syms_case_2
MOD_ORDER               =  1;
Demap_out_case_1 = demapper(rx_syms_case_1,mod_order_s,1);
Demap_out_case_2 = demapper(rx_syms_case_2,mod_order_s,1);
trel = poly2trellis(7, [171 133]); 
% viterbi decoder
rx_data_final_1= vitdec(Demap_out_case_1,trel,7,'trunc','hard');
rx_data_final_2 = vitdec(Demap_out_case_2,trel,7,'trunc','hard');

% rx_data is the final output corresponding to tx_data, which can be used
% to calculate BER

[number,ber3] = biterr(tx_interferer_data_a,rx_data_final_1);%for TXA data
[number,ber4] = biterr(tx_interferer_data_b,rx_data_final_2);%for TXB data

%reconstruction of interferer signal
[tx_data_a, reconstruct_A]=reconstruct_ofdm_tx(rx_data_final_1, pilots_A, mod_order_s, number_of_bits, N_SC, CP_LEN, SC_IND_DATA, SC_IND_PILOTS, N_OFDM_SYMS, trellis_end_length);
[tx_data_b, reconstruct_B]=reconstruct_ofdm_tx(rx_data_final_2, pilots_B, mod_order_s, number_of_bits, N_SC, CP_LEN, SC_IND_DATA, SC_IND_PILOTS, N_OFDM_SYMS, trellis_end_length);

%cancelling the interferer
H11_sv = repmat(h11_sv,N_OFDM_SYMS,1).';
H12_sv = repmat(h12_sv,N_OFDM_SYMS,1).';
H21_sv = repmat(h21_sv,N_OFDM_SYMS,1).';
H22_sv = repmat(h22_sv,N_OFDM_SYMS,1).';
% H11_sv = repmat(H1(1,1),N_OFDM_SYMS,64).';
% H12_sv = repmat(H1(1,2),N_OFDM_SYMS,64).';
% H21_sv = repmat(H1(2,1),N_OFDM_SYMS,64).';
% H22_sv = repmat(H1(2,2),N_OFDM_SYMS,64).';
canceller1 = H11_sv.*reconstruct_A+H12_sv.*reconstruct_B;
canceller2 = H21_sv.*reconstruct_A+H22_sv.*reconstruct_B;
syms_f_mat_pc_A = syms_f_mat_pc_A - canceller1;
syms_f_mat_pc_B = syms_f_mat_pc_B - canceller2;

syms_eq_mat_A = zeros(N_SC,N_OFDM_SYMS);
syms_eq_mat_B = zeros(N_SC,N_OFDM_SYMS);
for i = 1:N_OFDM_SYMS
    for j = 1:N_SC
         correction_fac = mul_vv(:,:,j)*[syms_f_mat_pc_A(j,i) syms_f_mat_pc_B(j,i)].';
         syms_eq_mat_A(j,i) = correction_fac(1);
         syms_eq_mat_B(j,i) = correction_fac(2);
    end
end

payload_syms_mat_A = syms_eq_mat_A(SC_IND_DATA, :);
payload_syms_mat_B = syms_eq_mat_B(SC_IND_DATA, :);



%% perform demodulate or demapping post combined symbols 
rx_syms_case_1 = reshape(payload_syms_mat_A,1,numel(payload_syms_mat_A));
rx_syms_case_2 = reshape(payload_syms_mat_B,1,numel(payload_syms_mat_B));
% plot the demodulated output rx_syms_case_1 and rx_syms_case_2
figure(4);subplot(2,2,3);
scatter(real(rx_syms_case_1), imag(rx_syms_case_1),'filled');
title(' Signal Space of received victim bits after hard SIC');
xlabel('I'); ylabel('Q');

figure(4);subplot(2,2,4);
scatter(real(rx_syms_case_2), imag(rx_syms_case_2),'filled');
title(' Signal Space of received victim bits after hard SIC');
xlabel('I'); ylabel('Q');


% FEC decoder for the rx_syms_case_1 and rx_syms_case_2
MOD_ORDER               =  1;
Demap_out_case_1 = demapper(rx_syms_case_1,MOD_ORDER,1);
Demap_out_case_2 = demapper(rx_syms_case_2,MOD_ORDER,1);
trel = poly2trellis(7, [171 133]); 
% viterbi decoder
rx_data_final_1= vitdec(Demap_out_case_1,trel,7,'trunc','hard');
rx_data_final_2 = vitdec(Demap_out_case_2,trel,7,'trunc','hard');

% rx_data is the final output corresponding to tx_data, which can be used
% to calculate BER

[number,ber5] = biterr(tx_victim_data_a,rx_data_final_1);%for TXA data
[number,ber6] = biterr(tx_victim_data_b,rx_data_final_2);%for TXB data

%plotting data
% snr_x = snr_start+[0:snr_iter-1];
% for x = 1:4
%     for y = 1:snr_iter
%         mean_mat_4x1(y,x) = mean(ber(:,y,x));
%         mean_mat_2x2(y,x) = mean(ber1(:,y,x));
%         var_mat_4x1(y,x) = var(ber(:,y,x));
%         var_mat_2x2(y,x) = var(ber1(:,y,x));
%     end
%     figure(6+2*(x-1));plot(snr_x,mean_mat_4x1(:,x));hold on;plot(snr_x,mean_mat_2x2(:,x));
%     title(sprintf(' Mean BER VS SNR for %s',modulations(x)));
%     xlabel('SNR(dB)');
%     ylabel('mean BER');
%     legend('BER for 4x1 system Q7.2','BER for 2x2 system Q7.3');
%     figure(6+2*(x)-1);plot(snr_x,var_mat_4x1(:,x));hold on;plot(snr_x,var_mat_2x2(:,x));
%     title(sprintf(' variance VS SNR for %s',modulations(x)));
%     xlabel('SNR(dB)');
%     ylabel('variance');
%     legend('variance for 4x1 system Q7.2','variance for 2x2 system Q7.3');
% end

%cancelling the interferer
% H11_sv = repmat(h11_sv,N_OFDM_SYMS,1).';
% H12_sv = repmat(h12_sv,N_OFDM_SYMS,1).';
% H21_sv = repmat(h21_sv,N_OFDM_SYMS,1).';
% H22_sv = repmat(h22_sv,N_OFDM_SYMS,1).';
% syms_f_mat_pc_A = syms_f_mat_pc_A - H11_sv.*syms_eq_mat_A-H12_sv.*syms_eq_mat_B;
% syms_f_mat_pc_B = syms_f_mat_pc_B - H21_sv.*syms_eq_mat_A-H12_sv.*syms_eq_mat_B;
% 
% syms_eq_mat_A = zeros(N_SC,N_OFDM_SYMS);
% syms_eq_mat_B = zeros(N_SC,N_OFDM_SYMS);
% for i = 1:N_OFDM_SYMS
%     for j = 1:N_SC
%          correction_fac = mul_vv(:,:,j)*[syms_f_mat_pc_A(j,i) syms_f_mat_pc_B(j,i)].';
%          syms_eq_mat_A(j,i) = correction_fac(1);
%          syms_eq_mat_B(j,i) = correction_fac(2);
%     end
% end