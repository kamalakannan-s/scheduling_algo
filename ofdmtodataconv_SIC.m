function [rx_syms_case] = ofdmtodataconv_SIC(raw_rx_dec,payload_ind,N_OFDM_SYMS,M)
% OFDM params
SC_IND_PILOTS           = [8 22 44 58];                           % Pilot subcarrier indices
SC_IND_DATA             = [2:7 9:21 23:27 39:43 45:57 59:64];     % Data subcarrier indices
N_SC                    = 64;                                     % Number of subcarriers
CP_LEN                  = 16;                                     % Cyclic prefix length
N_DATA_SYMS             = N_OFDM_SYMS * length(SC_IND_DATA);      % Number of data symbols (one per data-bearing subcarrier per OFDM symbol)

channel_coding = .5; % coding rate
trellis_end_length = 8; % bits for trellis to end

% Extract the payload samples (integral number of OFDM symbols following preamble)
for i = 1:M
rx_preamble_removed(i,:) = raw_rx_dec(i,payload_ind+(0:(CP_LEN+N_SC)*(N_OFDM_SYMS)-1));
payload_mat(:,:,i) = reshape(rx_preamble_removed(i,:),CP_LEN+N_SC,N_OFDM_SYMS);
% Remove the cyclic prefix
payload_mat_noCP(:,:,i) = payload_mat(CP_LEN+[1:N_SC], :,i);

end
syms_f_mat = fft(payload_mat_noCP, N_SC, 1);
for i = 1:M
    payload_syms_mat(:,:,i) = syms_f_mat(SC_IND_DATA, :,i);
    rx_syms_case(i,:) = reshape(payload_syms_mat(:,:,i),1,numel(payload_syms_mat(:,:,i)));
end

end
