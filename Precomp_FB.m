function [ basis, sample_points, runningtime] = Precomp_FB( n_r, R, c )
% Input:
%         n_r: number of points sampled in [0,c]
%         R: compact support radius
%         c: bandlimits
% Output:
%          basis.Phi_ns: bessel radial function J(R_{kq}\xi/c) for \xi in [0, c]
%          basis.ang_freqs: the associated angular frequencies (k in the paper)
%          basis.n_theta: the number of samples on concentric rings for polar Fourier transformation.
%          sample_points.r: positions in interval [0, c]
%          sample_points.w: weights
%          runningtime: time for generating the FB basis
%       Lingda Wang 04/21/2018

t1 = tic;

[ basis, sample_points ] = precomp_fb( n_r, R, c );

runningtime = toc(t1);