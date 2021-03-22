function raw_rx_dec_A=downsample_SIC(rx_vec_air_2A)
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
%% Decimate
if(INTERP_RATE == 1)
    raw_rx_dec_A = rx_vec_air_2A;
elseif(INTERP_RATE == 2)
    raw_rx_dec_A = filter(interp_filt2, 1, rx_vec_air_2A);
    raw_rx_dec_A = raw_rx_dec_A(1:2:end);
end


end
