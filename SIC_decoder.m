function [ber,rx_data_final] = SIC_decoder(H_sv,H_ss,H_vv,rank_s,raw_rx_dec,payload_ind,N_OFDM_SYMS,MOD_ORDER,precoder,M,tx_data_int)
K = 1.38e-23;
T = 300;
B = 20e6;
N0 = K*T*B;
[U_ss,S,V_ss] = svd(H_ss);
V = V_ss(:,1:rank_s);
G = H_sv*V_ss;
G_vv = H_vv*precoder;
%mul = G'*inv(N0*eye(M)+G_vv*G_vv'+G*G');
mul = inv(N0*eye(M)+G_vv'*G_vv+G'*G)*G';
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

%combining
for i = 1:N_OFDM_SYMS
    for j = 1:N_SC
        rx_mat =[];
        for k = 1:M
            rx_mat = [rx_mat syms_f_mat(j,i,k)];
        end
         correction_fac = mul*rx_mat.';
         syms_eq_mat(j,i,:) = correction_fac;
%          syms_eq_mat_B(j,i) = correction_fac(2);
    end
end
for i = 1:M
    payload_syms_mat(:,:,i) = syms_eq_mat(SC_IND_DATA, :,i);
    rx_data(i,:) = reshape(payload_syms_mat(:,:,i),1,numel(payload_syms_mat(:,:,i)));
end
trel = poly2trellis(7, [171 133]);
for i = 1:rank_s
figure(i);
scatter(real(rx_data(i,:)), imag(rx_data(i,:)),'filled');
title(' Signal Space of received victim bits after hard SIC');
xlabel('I'); ylabel('Q');
Demap_out_case(i,:) = demapper(rx_data(i,:),MOD_ORDER,1);
rx_data_final(i,:)= vitdec(Demap_out_case(i,:),trel,7,'trunc','hard');
[number,ber(i)] = biterr(tx_data_int(:,i).',rx_data_final(i,:));
end
end
