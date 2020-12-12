pi0 = 1.0;

tf = 5.0;
plow = -1.5;
pup = -0.5;

U = ChebopMarkovOneParameter(@Q,pi0,tf,plow,pup);

hold on
t=linspace(0,tf);
plot(t,U(t,-1));
plot(t,U(t,-0.5));


function out = Q(y)
out = y;
end