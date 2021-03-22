function [syms_precoded] = datatosymprecoded_SIC(rx_data_int, rank_int, MOD_ORDER, H_ss);
trel = poly2trellis(7, [171 133]);              % Define trellis
[U,S,V] = svd(H_ss);
precoder = V(:,1:rank_int);

for i = 1:rank_int
% Forward Error Correction
tx_code(i,:) = convenc(rx_data_int(i,:),trel);            % convultional encoder

% bits to signal space mapping these are you are x_k from the class
tx_syms(i,:) = mapping(tx_code(i,:), MOD_ORDER, 1);

end
syms_precoded = precoder*tx_syms;
end
