function [B] = B_matrix_creation_CU(B,L)
K = eye(2*L);               %twice the nodes because of UL and DL
K = [K zeros(2*L,1)];
[ dum,i] = size(K);
[ dum,j] = size(B);
c = [];
for col=1:i
    c =  [c repmat(K(:,col),1,j)];
end
d = repmat(B,1,i);
B = [d;c];
end