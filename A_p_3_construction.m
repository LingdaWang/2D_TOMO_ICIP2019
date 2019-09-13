function [A_p_3,b_p_3] = A_p_3_construction(Psi,x,y,index_coeff,k_max,Cov_E)

tmp1 = bsxfun(@times,Psi,x.');
tmp1 = mat2cell(tmp1,size(Psi,1),index_coeff);
tmp1 = cellfun(@(x) sum(x,2),tmp1,'UniformOutput', false);
tmp1 = cell2mat(tmp1);

tmp2 = bsxfun(@times,Psi',conj(y));
tmp2 = mat2cell(tmp2,index_coeff,size(Psi,1));
tmp2 = cellfun(@(x) sum(x,1),tmp2,'UniformOutput', false);
tmp2 = cell2mat(tmp2);
tmp4 = tmp1'*tmp1;

tmp5 = tmp2*tmp2';

%tmp3 = tmp1'*(Cov_E-tmp1*tmp2)*tmp2';
tmp3 = tmp1'*Cov_E*tmp2' - tmp4*tmp5;

%A_p_3 = zeros(2*k_max);

b_p_3 = matrix_shift(tmp3,k_max);

%{
for i=1:2*k_max
    A_p_3(:,i) = matrix_shift(tmp4*circshift(tmp5,[-i,0]),k_max);
end
%}
 A_p_3 = cell2mat(arrayfun(@(i) matrix_shift(tmp4*circshift(tmp5,[-i,0]),k_max),1:2*k_max,'un',0));












