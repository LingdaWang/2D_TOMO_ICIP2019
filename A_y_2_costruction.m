function [A_y_2,b_y_2] = A_y_2_costruction(x,C_P,Psi_Psi,Psi_b_Psi)

tmp = conj(x'.*C_P);
b_y_2 = sum(bsxfun(@times,tmp,Psi_b_Psi),2);
A_y_2 = bsxfun(@times,tmp*conj(Psi_Psi)*tmp',Psi_Psi);