function U = ChebopMarkovOneParameter(Q,pi0,tf,plow,pup)
%ChebopMarkovOneParameter approximate the solution bundle
%   consider the sulution bundle u(x,y) of the system of PDEs
%   du(x,y)/dx = Q(y)*u(x,y),
%   u(0,y) = pi0.
%   This system originates from the Kolmogorov Forward Equation (a system of ODEs)
%   dpi_p(t)/dt = Q_p*pi_p(t),
%   where the only independent variable is t, and p is a constant jet to be
%   defined, i.e., a parameter, that can be used to label the solutions.
%   If p is upgraded to the independent variable status, i.e., when we 
%   are looking for the solution bundle pi(t,p), renaming  
%   x:=t, y:=p and u:=u, we obtain the system of PDEs described at the
%   beginning. The found approximation of the solution bundle is called U.
%   Q is a function of a single variable, i.e., Q(y),
%   tf is such that the approximation U works for t=x in [0,tf],
%   pi0 must be a probability vector, i.e., 0<=pi0(i)<=1 and sum(pi0)=1,
%   plow and pup are such that plow<=p<=pup.

L = chebop2(@(x,y,u) diffx(u,1)-Q(y).*u, [0 tf plow pup]);
L.lbc = pi0;
rhs = zeros(length(pi0),1);
U = L\rhs;

end

