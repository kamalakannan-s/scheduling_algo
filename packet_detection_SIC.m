function [mimo_training_ind,payload_ind,rx_cfo_corr_t,fail]=packet_detection_SIC(raw_rx_dec_A,spacing,Total_cell,DO_APPLY_CFO_CORRECTION)
%% Correlate for LTS
LTS_CORR_THRESH=.8;
DO_APPLY_CFO_CORRECTION=0;
DO_APPLY_SFO_CORRECTION=0;
DO_APPLY_PHASE_ERR_CORRECTION=0;
fail = 0;
% LTS
lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
lts_t = ifft(lts_f, 64);

lts_corr = abs(conv(conj(fliplr(lts_t)), sign(raw_rx_dec_A)));
lts_corr = lts_corr(32:end-32);
lts_peaks = find(lts_corr > LTS_CORR_THRESH*max(lts_corr));
[LTS1, LTS2] = meshgrid(lts_peaks,lts_peaks);
[lts_last_peak_index,y] = find(LTS2-LTS1 == length(lts_t));
% Stop if no valid correlation peak was found
if(isempty(lts_last_peak_index))
    fprintf('No LTS Correlation Peaks Found!\n');
    fail = 1;
%     return;
end

% Set the sample indices of the payload symbols and preamble
% The "+32" here corresponds to the 32-sample cyclic prefix on the preamble LTS
% The "+192" corresponds to the length of the extra training symbols for MIMO channel estimation
mimo_training_ind = lts_peaks(max(lts_last_peak_index)) + 32;
mimo_training_ind = 662; %for trial
payload_ind = mimo_training_ind + Total_cell*spacing;

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
end

