function out = Q_10_12(y)
    mu = 1;
    out = [-2*y   2*y 0
             mu -mu-y y
              0    0 0];
end