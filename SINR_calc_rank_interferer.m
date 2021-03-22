function [SINR] = SINR_calc_rank_interferer(G_sv,G_vv,M,N0,m)

A = G_sv;
A(:,m) = [];
Z_sv = eye(M)*N0;
%for victim
Z_sv = Z_sv+G_vv*G_vv';
% % %         Z_sv = Z_sv + H(:,:,
% Z_sv = eye(M)*N0+H_vv*H_vv'+H_sv2*H_sv2';
W_sv = G_sv'*inv(G_sv*G_sv'+Z_sv);
num = abs(W_sv(m,:)*G_sv(:,m))^2;

den = W_sv(m,:)*(Z_sv+A*A')*W_sv(m,:)';
SINR = num/den;
end