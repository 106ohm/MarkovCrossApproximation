function pi = KolmogorovPoint2(Q, pi0, t, r, p, k)
%KOLMOGOROVPOINT 

% pi = expm(t * Q') * pi0;
% pi = expmv(t, Q', pi0);

% Make sure that pi0 is a row vector
pi0 = pi0(:)';

pi = k * pi0;

for j = 1 : length(r)
    fT = t*Q + p(j) * speye(size(Q));
    npi = r(j) * pi0 / fT;
    pi = pi - npi;
end

% Return a column vector
pi = pi(:);
pi = real(pi);

end

