function [A_x_2,b_x_2] = A_x_2_costruction(y,C_P,Psi_Psi,Psi_b_Psi)

tmp = conj(y'.*C_P);
b_x_2 = sum(bsxfun(@times,tmp,Psi_b_Psi),2);
A_x_2 = bsxfun(@times,tmp*conj(Psi_Psi)*tmp',Psi_Psi);

