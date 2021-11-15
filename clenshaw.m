function y = clenshaw(c, x)
% Clenshaw scheme for scalar-valued functions.
bk1 = 0*x; 
bk2 = bk1;
x = 2*x;
n = size(c,1)-1;
for k = (n+1):-2:3
    bk2 = c(k) + x.*bk1 - bk2;
    bk1 = c(k-1) + x.*bk2 - bk1;
end
if ( mod(n, 2) )
    tmp = bk1;
    bk1 = c(2) + x.*bk1 - bk2;
    bk2 = tmp;
end
y = c(1) + .5*x.*bk1 - bk2;
end
