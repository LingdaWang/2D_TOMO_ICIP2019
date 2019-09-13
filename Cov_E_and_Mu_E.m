function [Cov_E, Mu_E] = Cov_E_and_Mu_E(noisy_data,freqs,num_of_delta_alpha,n_r,Sigma)
FFT_noisy_data = zeros(2*n_r,size(noisy_data,2));
for i= 1:size(noisy_data,2)
    FFT_noisy_data(:,i) = nufft1(noisy_data(:,i),freqs);
end

FFT_noisy_data = reshape(FFT_noisy_data,2*n_r*(2*num_of_delta_alpha+1),[]);

Mu_E = mean(FFT_noisy_data,2);
Cov_E = FFT_noisy_data*FFT_noisy_data'./size(FFT_noisy_data,2);

for i=1:1:2*num_of_delta_alpha+1
    tmp3 = (i-1)*size(Sigma,1)+1;
    tmp4 = i*size(Sigma,1);
    Cov_E(tmp3:tmp4,tmp3:tmp4) = Cov_E(tmp3:tmp4,tmp3:tmp4) - Sigma;        
end


end