function [Mu_theoretical]=Mu_t_construction(Psi,a,C_P)
u_i = floor(size(C_P,1)/2);
Mu_theoretical = Psi*(a.*(C_P(u_i,:)'));
end