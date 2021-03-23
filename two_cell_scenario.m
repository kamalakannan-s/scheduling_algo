clear all;
clc;
close all;
c = 2;
L = 4;
M = 2;
slots = 500;
H_ref = zeros(M,M,c*L,c*L);
A = zeros(c*L,c*L);
A(4,5) = 1;
A(8,1) = 1;
% H_ref(:,:,4,5) = 1/sqrt(2)*complex(randn(M,M),randn(M,M));
% H_ref(:,:,8,1) = 1/sqrt(2)*complex(randn(M,M),randn(M,M));
N = 0.02*ones(1,c*L);
link_count_data = zeros(1,c*L);
link_count_fail = zeros(1,c*L);
count = 0;
for i = 1:2:c*L
    H_ref(:,:,i,i+1) = 1/sqrt(2)*complex(randn(M,M),randn(M,M));
    H_ref(:,:,i+1,i) = H_ref(:,:,i,i+1); %reciprocity
end


%finding the victim and it's interferer pair
v = [];
vs_pair = [];
for i = 1:c*L
    interference_array = A(i,:);
    interferer_num(i) = length(find(interference_array));
    if interferer_num(i) > 0
       v = [v i];
       vs_pair = [vs_pair find(interference_array)];
    end
end
V = length(v);
%defining channel between victim's rx and tx based on interferer's channel
%power
for i = 1:V
    H_ref(:,:,v(i),vs_pair(i)) = 10^(20/10)*H_ref(:,:,v(i)+mod(v(i),2)-1,v(i)+mod(v(i),2))+0.0002*complex(randn(M,M),randn(M,M));
%     H_ref(:,:,vs_pair(i),v(i))=  H_ref(:,:,v(i),vs_pair(i));
end
rate_mat = 1000*ones(c*L,V+1);
rank_mat = (M+1)*ones(c*L,V+1);
snr_min = 1000*ones(c*L,V+1);
H_vv = H_ref(:,:,v(1)+mod(v(1),2)-1,v(1)+mod(v(1),2));
s1 = 1:c*L;
    for i = 1:V
        s1(v(i)) = 0;
    end
s = find(s1);
S = length(s);
%finding standalone rank for links that are not victim
for i = 1:S
%     victim = 0;
    H(:,:,:) = H_ref(:,:,s(i),:);
    N0 = N(s(i)+(2*mod(s(i),2)-1));
    ind = s(i);
    interferer = find(A(ind,:));
    [standalone_rank(ind),rank_mat(ind,1),rate_mat(ind,1),snr_min(ind,1)] = standalone_rank_SIC_CU(H,N0,M,ind,interferer_num(ind),interferer);
end
%finding rank for interferer when victim is seleceted
%select rank for victim and use that rank for calculating 
%pre-SIC SINR of interferer
for i = 1:V
    H(:,:,:) = H_ref(:,:,v(i),:);
    N0 = N(v(i)+(2*mod(v(i),2)-1));
    ind = v(i);
    interferer = find(A(ind,:));
    [standalone_rank(ind),rank_mat(ind,1),rate_mat(ind,1),snr_min(ind,1)] = standalone_rank_SIC_CU(H,N0,M,ind,interferer_num(ind),interferer);
    k = length(interferer);
    for j = 1:k
        current = interferer(j);
        [U,S,prec] = svd(H(:,:,ind+(2*mod(ind,2)-1)));
        G_vv = H(:,:,ind+(2*mod(ind,2)-1))*prec(:,1:standalone_rank(ind));
        H_sv = H(:,:,current);
        H_ss = H_ref(:,:,current+mod(current,2)-1,current+mod(current,2));
        [standalone_rank(current),rank_mat(current,i+1),rate_mat(current,i+1),snr_min(current,i+1)] = standalone_rank_SIC_CU_interferer(H_ss,H_sv,G_vv,N0,M);
%     standalone_rank(ind) = standalone_rank_SIC_CU_victim(H,N0,M,ind,interferer_num(ind),interferer,H_ref,standalone_rank);
    end

end

for i = 1:V
rate_mat(:,i+1) = min(rate_mat(:,1),rate_mat(:,i+1));
rank_mat(:,i+1) = min(rank_mat(:,1),rank_mat(:,i+1));
end

%assuming 2 users per cell and 2 cells in total
c = 2;
l = 2;
B = 1;
B = B_matrix_creation_CU(B,l);
B = B(2:end,:);
for i = 1:c-1
    B = B_matrix_creation_CU(B,l);
end

iter = length(B); %always return the largest dimension of matrix, column is always large
w = 130*ones(c*L,1);
for slots_ind = 1:slots

% Calculate weights
% % % % r = poissrnd(0.0446, 2*L, 1); % (6 * 10^6 / 2000)/(8400 * 8)
% % % % w = min(r, 1);
% % % % if slots_ind > 1
% % % %     w = w + w_new;
% % % % end
% % % % w_new = zeros(2*L, 1);
    
for i = 1:iter
    b = B(:,i);
    rate_sel = rate_mat(:,1); % assign standalone rate in the beginning
    selected_ind = find(b');
    victim_sel = ismember(v,selected_ind);
    for j = 1:V %checking whether any of the victim is selected or not
        if victim_sel(j)
            rate_sel = min(rate_sel,rate_mat(:,j)); %jth victim is selected in the scheduling combination under consideration
        end
%         rank_mat(:,j+1) = min(rank_mat(:,j),rank_mat(:,1));
%         snr_min(:,j+1) = min(snr_min(:,j+1),snr_min(:,1));
    end
    obj(i) = sum(b.*w.*rate_sel);
end
  [dummy B_ind] =  max(obj);
  b = B(:,B_ind);
%selecting rank for the scheduled links
rank_sel = rank_mat(:,1); % assign standalone rate in the beginning
selected_ind = find(b');
victim_sel = ismember(v,selected_ind);
    for j = 1:V %checking whether any of the victim is selected or not
        if victim_sel(j)
            rank_sel = min(rank_sel,rank_mat(:,j)); %jth victim is selected in the scheduling combination under consideration
        end
    end
    rank_sel = rank_sel.*b;
%rank_sel is the rank of the links and b is the scheduled link for the slot
%call the cells(both tx and rx immediately) to check for successful
%transmission
%if successful, update the corresponding weigth
%now go to next slot scheduling
link = find(b');
cell = length(link);
for i = 1:cell
    cell_number = ceil(link(i)/l);
    [lia,loc] = ismember(link(i),v);
    rank = rank_sel(link(i));
    mod_order = 1;
    H_vv = H_ref(:,:,link(i)+mod(link(i),2)-1,link(i)+mod(link(i),2)); %check once more
    if lia %checking whether the scheduled link is victim or not
        fprintf("victim is selected\n");
        if ismember(vs_pair(loc),link)
            fprintf("both the interferer and victim are selected\n");
            rank_int = rank_sel(vs_pair(loc));
            H_sv = H_ref(:,:,link(i),vs_pair(loc));
            H_ss = H_ref(:,:,vs_pair(loc)+mod(vs_pair(loc),2)-1,vs_pair(loc)+mod(vs_pair(loc),2));
            cell_number_int = ceil(vs_pair(loc)/l);
            mod_order_int = 1;
            [success(i),N(link(i))] = SIC_PHY(rank,rank_int,mod_order,mod_order_int,cell_number,cell_number_int,H_ss,H_vv,H_sv);
        else
            [success(i),N(link(i))] = SVD_PHY(rank,mod_order,cell_number,H_vv);
        end
    else
        [success(i),N(link(i))] = SVD_PHY(rank,mod_order,cell_number,H_vv);
    end
    if success
        %update weight of link(i) position
        w_new(link(i)) = 0;
        link_count_data(link(i)) = link_count_data(link(i))+1;
        w(link(i)) = w(link(i))-1;
    else
        w_new(link(i)) = w(link(i)) + 1;
        count = count+1;
        link_count_fail(link(i)) = link_count_fail(link(i))+1;
    end
end
end
figure(5);plot(link_count_data);
xlabel("Links");
ylabel("packets transmitted successfully");
title("packets transmitted across links for 500 slots");
figure(6);plot(link_count_fail);
xlabel("Links");
ylabel("packets dropped");
title("packets dropped across links for 500 slots");
figure(7);plot(abs(rate_mat(:,1)));
xlabel("Links");
ylabel("Standalone rate of links");
title("Standalone rate of links for realization");