function Sigma = Sigma_constr(noise_var,R,n_r,freqs)
tmp_test = diag(noise_var.*ones(1,2*R));
tmp_test2 = [];
tmp_test3 = [];
for i = 1:2*R
    tmp_test2 = [tmp_test2,nufft1(tmp_test(:,i),freqs)];
end
for i = 1:2*n_r
    tmp_test3 = [tmp_test3,nufft1(tmp_test2(i,:)',freqs)];
end
Sigma = tmp_test3';
end