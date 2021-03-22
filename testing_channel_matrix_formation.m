clear all
close all
clc
% for i = 1:4
%     a(i,:) = i*ones(1,10);
% end
% M = 2;
% for i = 1:4
%     h((floor((i-1)/M))+1,mod((i-1),M)+1,:) = a(i,:);
% end
% for i = 1:4
% b(:,:,i) = reshape(a(i,:),2,5);
% end
B = 1;
% B = [B zeros(2,1)];
B = B_matrix_creation_CU(B,2);
% B = B_matrix_creation_CU(B,4);
% % % M = 2;
% % % for i = 1:5
% % %     H(:,:,i) = i*[1 2;3 4];
% % % end
% % % H_m = mean(H,3);