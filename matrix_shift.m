function [result] = matrix_shift(A,k_max)

% tmp1 = mat2cell(A,size(A,1),ones(1,size(A,1)));
% tmp2 = cellfun(@(x) circshift(x,))
result = sum(cell2mat(arrayfun(@(i) circshift(A(:,i),1-i),1:2*k_max+1,'UniformOutput',0)),2);
result = flipud(result(2:end));
