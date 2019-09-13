function [A_p_1, tmp_x_p_1] = A_p_1_costruction(Psi,x,index_coeff,k_max)

tmp1 = bsxfun(@times, Psi, x.');
tmp2 = fliplr(mat2cell(tmp1,size(Psi,1),index_coeff));
tmp3 = cell2mat(cellfun(@(x) sum(x,2), tmp2, 'UniformOutput', false));
A_p_1 = [tmp3(:,k_max+2:end),tmp3(:,1:k_max)];
tmp_x_p_1 = tmp3(:,k_max+1);
