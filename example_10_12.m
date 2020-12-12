pi0 = [1.0;0;0];

tf = 5.0;
plow = 0.5;
pup = 1.5;

U = ChebopMarkovOneParameter(@Q,pi0,tf,plow,pup);

function out = Q(y)
mu = 1;
out = [-2*y   2*y 0
         mu -mu-y y
          0    0 0];
end