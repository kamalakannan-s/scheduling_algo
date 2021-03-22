function [SINR] = SINR_calc_rank(H_vv,H,m,M,N0,s,interferer);

A = H_vv;
A(:,m) = [];
Z_sv = eye(M)*N0;
%for victim
% % % if s > 0
% % %     count = length(interferer);
% % %     for i = 1:s
% % %         Z_sv = Z_sv + H(:,:,
% Z_sv = eye(M)*N0+H_vv*H_vv'+H_sv2*H_sv2';
W_sv = H_vv'*inv(H_vv*H_vv'+Z_sv);
num = abs(W_sv(m,:)*H_vv(:,m))^2;

den = W_sv(m,:)*(Z_sv+A*A')*W_sv(m,:)';
SINR = num/den;
end