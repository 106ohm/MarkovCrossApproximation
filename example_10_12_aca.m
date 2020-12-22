pi0 = [1.0;0;0];

tf = 5.0;
plow = 0.5;
pup = 1.5;

t = linspace(0, tf, 100);

% pi is a matrix whose rows are pi(t_i) for t_i in [0,tf]
% pi = KolmogorovODE(Q(lambda), pi0, t);

lambda = linspace(0.5, 1.5, 100);

%U = ChebopMarkovOneParameter(@Q,pi0,tf,plow,pup);
n = [ length(t), length(lambda) ];
Afiber = @(j,i) example_10_12_fiber(j, i, t, lambda, @Q, pi0);
U = aca_nd(n, Afiber, 1e-6);

function out = Q(y)
    mu = 1;
    out = [-2*y   2*y 0
             mu -mu-y y
              0    0 0];
end