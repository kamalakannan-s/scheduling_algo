function SINR = SINR_calc(H_sv1,H_sv2,H_vv,m,M,N0)

Z_sv = eye(M)*N0+H_vv*H_vv'+H_sv2*H_sv2';
W_sv = H_sv1'*inv(H_sv1*H_sv1'+Z_sv);
num = abs(W_sv(m,:)*H_sv1(:,m))^2;
A = H_sv1;
A(:,m) = [];
den = W_sv(m,:)*(Z_sv+A*A')*W_sv(m,:)';
SINR = num/den;
end