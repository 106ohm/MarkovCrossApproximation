function pi = KolmogorovPoint(Q, pi0, t)
%KOLMOGOROVPOINT 

pi = expm(t * Q') * pi0;

end

