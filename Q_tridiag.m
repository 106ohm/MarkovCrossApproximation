function Q = Q_tridiag(y, m)
%Q_TRIDIAG 

u = y(1) * (m : -1 : 1)'; u(1) = 0;
l = y(2) * ones(m-1, 1);  l(m) = 0;

d = zeros(m,1);
d(1) = -u(2);
d(m) = -l(m-1);
d(2:m-1) = -(l(1:m-2) + u(3:m));

d(m) = 0;
l(m-1) = 0;

Q = spdiags([ l, d, u ], -1:1, m, m);

end

