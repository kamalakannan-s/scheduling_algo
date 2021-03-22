function [opt_rank] = standalone_rank_SIC_CU_victim(H,N0,M,ind,s,interferer,H_ref,standalone_rank)
H_vv = H(:,:,ind+(2*mod(ind,2)-1));
[U,S,V] = svd(H_vv);
for j = 1:M
for i = 1:j
    G_vv = H_vv*V(:,1:i);
    SINR(i) = SINR_calc_rank_victim(G_vv,H,i,M,N0,s,interferer,H_ref,standalone_rank);
end
SINR_db = 10*log10(SINR)
SINR_sum(j) = sum(SINR_db);
end
[max_val opt_rank] = max(SINR_sum);
end