function [ Psi, runningtime] = Psi_Construction( basis, num_of_delta_alpha, delta_alpha,k_max)

% Description: Construct the Psi martix
% Input:
%   basis: FB basis
%   alpha: angle range we cover for each sample theta
%   smallest_interval: delta_alpha we choose
%   k_max: highest angular frequency
% Output:
%   Psi: Psi matrix
%   runntingtime: time for running the function
%   
%   Lingda Wang 04/21/2018


t1 = tic;
%% construct a line of Psi matrix with out phase

Phi_flip = fliplr(reshape(basis.Phi_ns,1,[]));

Phi = cat(2,Phi_flip(1:end-1),reshape(basis.Phi_ns,1,[]));

%num_of_delta_alpha = ceil(alpha./smallest_interval);

k = -k_max:1:-1;

temp1 = (-1).^k;

temp1 = [temp1,ones(1,1+k_max)];

Phi = cellfun(@times,num2cell(temp1),Phi,'UniformOutput', false);
Phi_2 = cellfun(@(x) flipud(x),Phi,'UniformOutput', false);

Psi = [];

%% construct the whole Psi matrix based on Phi and also adding the phase

for i = -num_of_delta_alpha:1:num_of_delta_alpha
    
    %delta_alpha_1 = exp(sqrt(-1).*(-k_max:1:k_max).*i.*delta_alpha);
    delta_alpha_1 = exp(sqrt(-1).*(-k_max:1:k_max).*i.*delta_alpha);
    %delta_alpha_1 = exp(sqrt(-1).*(-k_max:1:k_max).*i.*0);
   % delta_alpha_2 = exp(sqrt(-1).*(-k_max:1:k_max).*(i.*delta_alpha+pi));
    delta_alpha_2 = exp(sqrt(-1).*(-k_max:1:k_max).*(i.*delta_alpha+pi));
    temp = cellfun(@times,Phi,num2cell(delta_alpha_1) ,'UniformOutput', false);
    temp2 = cellfun(@times,Phi_2,num2cell(delta_alpha_2) ,'UniformOutput', false);
    Psi = cat(1, Psi, temp2);
    Psi = cat(1, Psi, temp);
end

Psi = cell2mat(Psi);

runningtime = toc(t1);
