%% Model fixed parameters
ncycles = 7;
mu= 0.5; 
lambda = 1e-4;
D85 = 2;
D15 = 0.1;

% approximate deterministic distribution with CDF delta(t-D85)
% as an Erlang with n_erlang_d85 states,
% and similarly for  the CDF delta(t-D15)
n_erlang_d85 = 5;
n_erlang_d15 = 3;

nstates = (ncycles-1)*(1+n_erlang_d85+n_erlang_d15+1)+1+1;

pi0 = zeros(nstates,1);
pi0(1) = 1;

Q = @(y1,y2)infgen(nstates, ncycles, mu, lambda, y1, y2, n_erlang_d85, n_erlang_d15);

%% Model variable parameters (time and D15)
tf = 5.0;
plow = 0.1;
pup = 0.5;

t = linspace(0, tf, 100);

D15 = plow;
D85 = plow;
% pi is a matrix whose rows are pi(t_i) for t_i in [0,tf]
pi = KolmogorovODE(Q(D15, D85), pi0, t);

D15 = linspace(plow, pup, 100);
D85 = linspace(plow, pup, 100);

%U = ChebopMarkovOneParameter(@Q,pi0,tf,plow,pup);


%% Measure of interest
r = [ones(nstates-ncycles,1);zeros(ncycles,1)];% Reliability at time t

%% ACA related definitions
Afiber = @(j,i) dot(r, KolmogorovODE(  Q( D15(i(2)), D85(i(3)) )  ), pi0, t ); % TODO fix this

n = [ length(t), length(D15), length(D85) ];

U = aca_nd(n, Afiber, 1e-6);

%% Infinitesimal Generator Matrix
function Q = infgen(nstates, ncycles, mu, lambda, y1, y2, n_erlang_d85, n_erlang_d15)
% the last ncycles states are (absorbing) failure states

% the states of each cycle (for cycle<ncycles) are indexed as follows:
% 1 is C>15,
% from 2 to 1+n_erlang_d85 are D>15,
% from 2+n_erlang85 to 1+n_erlang_d85+n_erlang_d15 are D<15,
% nstates-ncycles+cycle is F.
% the state of the last cycle are indexes as follows:
% 1 is C>15, 2 is F

R = sparse(0,0);
for cycle = 1 : ncycles-1
    base_index = (cycle-1)*(1+n_erlang_d85+n_erlang_d15+1);
    Rcycle = sparse(1+n_erlang_d85+n_erlang_d15+1,nstates);
    Rcycle(1,base_index+2) = lambda;
    Rcycle(2,base_index+1) = mu;
    for i = 1 : n_erlang_d85
        Rcycle(1+i,base_index+1+i+1) = n_erlang_d85/y2;
        Rcycle(1+i,base_index+1) = mu;
    end
    for i = 0 : n_erlang_d15-1
        Rcycle(1+n_erlang_d85+i,base_index+1+n_erlang_d85+i+1) = n_erlang_d15/y1;
        Rcycle(1+n_erlang_d85+i,base_index+1+n_erlang_d85+n_erlang_d15+1) = mu;
    end
    Rcycle(1+n_erlang_d85+n_erlang_d15,nstates-ncycles+cycle) = n_erlang_d15/y1;
    R = [R; Rcycle];
end
R(end+1,end) = lambda;
R(end+1,:) = 0;

Q = (R-diag(R*ones(nstates,1)));
end