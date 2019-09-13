function [ coeff, runningtime ] = Coeff_FB(P, R, noise_variance, basis, sample_points, num_pool)

%Description of £ºCoeff_FB

% Computes Fourier-Bessel expansion coefficients 'coeff'
%Input:
%	P: image.
%	c: band limit
%	R: compact support radius
%	noise_variance: estimated noise variance
%Output
%	coeff: Fourier-Bessel expansion coefficients of input image
%   Lingda Wang 04/21/2018

t1 = tic;

[ ~, coeff, ~, ~, ~, ~ ] = jobscript_FFBsPCA( P, R, noise_variance, basis, sample_points, num_pool );

runningtime = toc(t1);
