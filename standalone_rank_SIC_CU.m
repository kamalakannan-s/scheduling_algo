function [opt_rank rank_mat rate_val SINR_min] = standalone_rank_SIC_CU(H,N0,M,ind,s,interferer)
H_vv = H(:,:,ind+(2*mod(ind,2)-1));
[U,S,V] = svd(H_vv);
for j = 1:M
    clear SINR_db;
for i = 1:j
    G_vv = H_vv*V(:,1:i);
    SINR(i) = SINR_calc_rank(G_vv,H,i,M,N0,s,interferer);
end
SINR_db = 10*log10(SINR);
SINR_sum(j) = sum(SINR_db)-(j-1)*6;
end
[max_val opt_rank] = max(SINR_sum);
rank_mat = opt_rank;
G_vv = H_vv*V(:,1:opt_rank);
for i = 1:opt_rank
    SINR_opt(i) = SINR_calc_rank(G_vv,H,i,M,N0,s,interferer);
end
rate = 0.8*log2(1+SINR_opt/1.333);
rate_val = sum(rate);
SINR_db_opt = 10*log10(SINR_opt);
SINR_min = min(SINR_db_opt);
end