function [ a, runningtime ] = a_Construction(coeff)
%Description: construct the a vector
%Input: 
%   coeff: Fourier-Bessel expansion coefficients of input image
%Output:
%   a: vector of full coeff in the order of [a_{-k_{max}q};...;a_{k_{max}q}]
%   runntingtime: total time of running
%   Lingda Wang 04/21/2018

t1 = tic;

a1 = flipud(coeff); % flip the coeff updown

a1 = conj(cell2mat(a1)); % a_{-k,q} = conj(a_{kq})

p0 = size(coeff{1}, 1); 

a_minus = a1(1:end-p0); % coeff for the part of k<0

a = [ a_minus; cell2mat(coeff) ]; % full coeff 

runningtime = toc(t1);