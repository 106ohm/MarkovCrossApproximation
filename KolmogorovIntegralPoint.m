function pi = KolmogorovIntegralPoint(Q, pi0, t, r, p, k)
%KOLMOGOROVPOINT 

% Make sure that pi0 is a row vector
pi0 = pi0(:)';

pi = k * pi0;

for j = 1 : length(r)
    fT = t*Q + p(j) * speye(size(Q));
    npi = - ( r(j) / p(j) ) * pi0 / fT;
    pi = pi - npi;
end

% Return a column vector
pi = t * pi(:);
pi = real(pi);

end

