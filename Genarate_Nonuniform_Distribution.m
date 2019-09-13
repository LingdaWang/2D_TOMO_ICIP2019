function [ pmf_theta, num_of_sample,rand_num] = Genarate_Nonuniform_Distribution( k_max, num_of_sample )
%% Generate the non_uniform pmf
rng(1);
tmp1 = 2*k_max+1;% num of pmf's
rand_num = ceil(rand(tmp1,1)*tmp1);% make sure each pmf is not travial.
%rand_num = ones(tmp1,1).*ceil(num_of_sample/tmp1);
pmf_theta = rand_num./sum(rand_num,1);% compute the pmf
rand_num = ceil(pmf_theta.*num_of_sample);% count how many samples each angles have
num_of_sample = sum(rand_num);% re-count the num of samples
pmf_theta = rand_num./num_of_sample; % re-compute the pmf
