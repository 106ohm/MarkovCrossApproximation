pi0 = [1.0;0;0];

tf = 5.0;
plow = 0.5;
pup = 1.5;

t = linspace(0, tf, 100);

lambda = 0.1;

% pi is a matrix whose rows are pi(t_i) for t_i in [0,tf]
pi = KolmogorovODE(Q(lambda), pi0, t);

lambda = linspace(0.5, 1.5, 100);

%U = ChebopMarkovOneParameter(@Q,pi0,tf,plow,pup);

e1 = [ 1 ; 0 ; 0 ];

Aelem = @(i,j) dot(e1, KolmogorovPoint(Q(lambda(j)), pi0, t(i)));
Arow  = @(i) arrayfun(@(j) Aelem(i, j), 1 : length(lambda))';
Acol  = @(j) KolmogorovODE(Q(lambda(j)), pi0, t) * e1;

[U, V] = aca_2d(length(t), length(lambda), Aelem, Arow, Acol);

function out = Q(y)
    mu = 1;
    out = [-2*y   2*y 0
             mu -mu-y y
              0    0 0];
end