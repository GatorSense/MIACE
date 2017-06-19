function [ace_data,mu,siginv] = ace_det(hsi_data,tgt_sig,mu,siginv,targflag)
% input anything for targflag if target signature is pulled from data


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

ace_data = A./(B.*C);

end