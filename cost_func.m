function [cost_func_value,cost1,cost2,cost3,cost4] = cost_func(lambda_1,lambda_2,rho,u,x,y,C_P,Mu_E,Cov_E,Psi,u_i)

cost1 = lambda_1/2*norm(Psi*(x.*(C_P(u_i,:)'))-Mu_E)^2;
cost2 = lambda_1/2*norm(Psi*(y.*(C_P(u_i,:)'))-Mu_E)^2;
cost3 = lambda_2/2*norm(Psi*(x*y'.*C_P)*Psi'-Cov_E,'fro')^2;
cost4 = rho/2*norm(x-y+u)^2;

cost_func_value = cost1+cost2+cost3+cost4;
