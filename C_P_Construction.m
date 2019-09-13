function [ index_C_P, C_P, runningtime ] = C_P_Construction( FFT_pmf,coeff_structure )
%Description: Construct the C_P matrix which is consisted of FFT(pmf_theta)
%Input: 
%   coeff_structure: Fourier-Bessel expansion coefficients structure
%   FFT_pmf: fft(the probability mass function of projections taken at theta)
%   k_max:highest angular frequency
%Output:
%   C_P:C_P matrix
%   index: store the information of which FFT(pmf_theta) that C_P_{ij} refers to
%   runntingtime: total time of running
%   Lingda Wang 04/21/2018

 t1 = tic;
 
 k_max = length(coeff_structure)-1;
 %% Compute the matrix index_C_P and index_C_P_{ij} = k, where C_P_{ij} = FFT(pmf_theta)_k 

 k_minus = flipud(coeff_structure);
 
 k_vec = cat(1,k_minus(1:end-1),coeff_structure);
 
 k_vec = cellfun(@times,k_vec,num2cell((-k_max:1:k_max)'),'UniformOutput',false);
 
 k_vec =int64(round(cell2mat(k_vec)));
 
 index_C_P = -(k_vec-k_vec');
 %% pad FFT(pmf_theta) with 0 if size(FFT_pmf,1) < 2*k_max + 1
 
 if size(FFT_pmf,1) < 2*k_max + 1

     FFT_pmf = [FFT_pmf;zeros(2*k_max + 1-size(FFT_pmf,1),1)];
     
 end
 
 %% give FFT(pmf_theta)_k, where k<0
 
 FFT_pmf_minus = conj(flipud(FFT_pmf));
 FFT_pmf = [FFT_pmf_minus(end-2*k_max:end-1);FFT_pmf(1:2*k_max+1)];
 
 %% Get the C_P matrix
 
 C_P = FFT_pmf( index_C_P+2*k_max+1 );
 
 runningtime = toc(t1);