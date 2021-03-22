function [preamble,spacing] = preamble_generator(M,cell,Total_cell)

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

preamble_legacy(1,:) = [sts_t_rep, lts_t(33:64), lts_t, lts_t];
preamble_legacy(2,:) = [circshift(sts_t_rep, [0, TX_SPATIAL_STREAM_SHIFT]), zeros(1, 160)];


lts_t_data = [lts_t(33:64), lts_t];
for i = 1:M
    preamble_mimo(i,:) = [zeros(1,96*(i-1)),lts_t_data,zeros(1,96*(M-i))];
    preamble(i,:) = [preamble_legacy(i,:),zeros(1,length(preamble_mimo(i,:))*(cell-1)),preamble_mimo(i,:),zeros(1,length(preamble_mimo(i,:))*(Total_cell-cell))];
end
spacing = length(preamble_mimo(1,:));
end