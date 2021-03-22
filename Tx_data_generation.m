function [tx_data,tx_code,preamble_legacy_A,preamble_legacy_B,preamble_mimo_A,preamble_mimo_B] = Tx_data_generation(M,rank,sym_size,Mod_order)
N_OFDM_SYMS             = sym_size;         % Number of OFDM symbols
MOD_ORDER               =  Mod_order;          % Modulation order in power of 2 (1/2/4/6 = BSPK/QPSK/16-QAM/64-QAM)
TX_SCALE                = 1.0;         % Scale for Tx waveform ([0:1])
% OFDM params
SC_IND_PILOTS           = [8 22 44 58];                           % Pilot subcarrier indices
SC_IND_DATA             = [2:7 9:21 23:27 39:43 45:57 59:64];     % Data subcarrier indices
N_SC                    = 64;                                     % Number of subcarriers
CP_LEN                  = 16;                                     % Cyclic prefix length
N_DATA_SYMS             = N_OFDM_SYMS * length(SC_IND_DATA);      % Number of data symbols (one per data-bearing subcarrier per OFDM symbol)
channel_coding = .5; % coding rate
trellis_end_length = 8; % bits for trellis to end
%to be removed
TX_SPATIAL_STREAM_SHIFT=3;
% STS
sts_f = zeros(1,64);
sts_f(1:27) = [0 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0];
sts_f(39:64) = [0 0 1+1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0];
sts_t = ifft(sqrt(13/6).*sts_f, 64);
sts_t = sts_t(1:16);

% LTS
lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
lts_t = ifft(lts_f, 64);

sts_t_rep = repmat(sts_t, 1, 30);

preamble_legacy_A = [sts_t_rep, lts_t(33:64), lts_t, lts_t];
preamble_legacy_B = [circshift(sts_t_rep, [0, TX_SPATIAL_STREAM_SHIFT]), zeros(1, 160)];

preamble_mimo_A = [lts_t(33:64), lts_t, zeros(1,96)];
preamble_mimo_B = [zeros(1,96), lts_t(33:64), lts_t];
%upto this
pilots_A= [1 1 -1 1].';
pilots_B= [0 0 0 0].';

%% Generate a payload of random integers
number_of_bits= (N_DATA_SYMS * MOD_ORDER - 2*trellis_end_length) * channel_coding;

    tx_data = randi(2, 1, number_of_bits) - 1;
    
    tx_data = double([tx_data zeros(1,trellis_end_length) ]);    % 8 bits padding
    trel = poly2trellis(7, [171 133]);              % Define trellis
    tx_code = convenc(tx_data,trel);            % convultional encoder

end
