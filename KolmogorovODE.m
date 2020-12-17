function pi = KolmogorovODE(Q, pi0, t)
%KOLMOGOROVODE solve the Forward Kolmogorov system of ODEs
%   Q is the infinitesimal generator matric
%   pi0 is the initial condition, i.e., a probability vector
%   the equations is integrate from 0 to tf

[t,pi] = ode45(@(t,pi) Q'*pi, t, pi0);

end

