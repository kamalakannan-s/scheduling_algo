function [ifft_in_mat] = payload_mat_gen_SIC(syms_precoded,n,N_OFDM_SYMS)

%OFDM data
pilots(:,1)= [1 1 -1 1].';
pilots(:,2)= [0 0 0 0].';
SC_IND_PILOTS           = [8 22 44 58];                           % Pilot subcarrier indices
SC_IND_DATA             = [2:7 9:21 23:27 39:43 45:57 59:64];     % Data subcarrier indices
N_SC                    = 64;                                     % Number of subcarriers
CP_LEN                  = 16;

% Reshape the symbol vector to a matrix with one column per OFDM symbol,
tx_syms_mat = reshape(syms_precoded, length(SC_IND_DATA), N_OFDM_SYMS);


% Repeat the pilots across all OFDM symbols
pilots_mat = repmat(pilots(:,n), 1, N_OFDM_SYMS);


% Construct the IFFT input matrix
ifft_in_mat = zeros(N_SC, N_OFDM_SYMS);

% Insert the data and pilot values; other subcarriers will remain at 0
ifft_in_mat(SC_IND_DATA, :)   = tx_syms_mat;
ifft_in_mat(SC_IND_PILOTS, :) = pilots_mat;
end