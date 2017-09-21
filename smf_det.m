function [smf_data,mu,siginv] = smf_det(hsi_data,tgt_sig,mu,siginv,target_flag)
% SMF Spectral Matched Filter
%
% Syntax:  [smf_data,mu,siginv] = smf_det(hsi_data,tgt_sig,mu,siginv)
%
% Inputs:
%   hsi_data - DxN array of N data points of dimensionality D
%   tgt_sig - Dx1 vector containing the target signature
%   mu - Dx1 vector containing the background mean vector, if empty,
%       computed as mean of all hsi_data
%   siginv - DxD matrix containing the inverse background covariance, if
%       empty, computed from all hsi_data
%   target_flag - flag indicating whether mean should be subtracted from
%   	target signatures or not, set to anything if mean should be subtracted.
% Outputs:
%   smf_data - Nx1 vector of SMF confidence values corresponding to each
%       test point
%   mu - Dx1 vector containing the background mean vector
%   siginv - DxD matrix containing the inverse background covariance
%
% University of Florida, Electrical and Computer Engineering
% Email Address: azare@ufl.edu
% Latest Revision: September 21, 2017
% This product is Copyright (c) 2017 University of Florida
% All rights reserved.



if isempty(mu)
    mu = mean(hsi_data,2);
end
if isempty(siginv)
    siginv = pinv(cov(hsi_data'));
end

if (nargin <5)
s = tgt_sig;
else
s = tgt_sig - mu;
end

z = bsxfun(@minus,hsi_data,mu);

st_siginv = s'*siginv;
st_siginv_s = s'*siginv*s;


A = sum(st_siginv*z,1);
B = sqrt(st_siginv_s);
C = sqrt(sum(z.*(siginv*z),1));

smf_data = A./(B);

end