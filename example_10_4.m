pi0 = [1.0;0;0;0];

tf = 5.0;
plow = 0.5;
pup = 1.5;

lambda = [0.75, 1.25];
% pi is a matrix whose rows are pi(t_i) for t_i in [0,tf]
pi = KolmogorovODE(Q(lambda),pi0,tf);

%U = ChebopMarkovOneParameter(@Q,pi0,tf,plow,pup);

function out = Q(y)
out = [-2.5     1  1.5    0
          0 -y(2)    0  y(2)
          0     0 -y(1) y(1)
          0     0    0    0];
end