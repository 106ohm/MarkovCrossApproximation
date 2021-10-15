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

Q = @(y)infgen(nstates, ncycles, mu, lambda, D85, y, n_erlang_d85, n_erlang_d15);

%% Model variable parameters (time and D15)
tf = 5.0;
plow = 0.1;
pup = 0.5;

t = linspace(0, tf, 100);

D15 = 0.1;
% pi is a matrix whose rows are pi(t_i) for t_i in [0,tf]
pi = KolmogorovODE(Q(D15), pi0, t);

D15 = linspace(plow, pup, 100);

%U = ChebopMarkovOneParameter(@Q,pi0,tf,plow,pup);


%% Measure of interest
r = [ones(n-ncycles,1);zeros(ncycles,1)];% Reliability at time t

%% ACA related definitions
Aelem = @(i,j) dot(r, KolmogorovPoint(Q(D15(j)), pi0, t(i)));
Arow  = @(i) arrayfun(@(j) Aelem(i, j), 1 : length(D15))';
Acol  = @(j) KolmogorovODE(Q(D15(j)), pi0, t) * r;

[U, V] = aca_2d(length(t), length(D15), Aelem, Arow, Acol);

%% Infinitesimal Generator Matrix
function Q = infgen(nstates, ncycles, mu, lambda, D85, y, n_erlang_d85, n_erlang_d15)
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
        Rcycle(1+i,base_index+1+i+1) = n_erlang_d85/D85;
        Rcycle(1+i,base_index+1) = mu;
    end
    for i = 0 : n_erlang_d15-1
        Rcycle(1+n_erlang_d85+i,base_index+1+n_erlang_d85+i+1) = n_erlang_d15/y;
        Rcycle(1+n_erlang_d85+i,base_index+1+n_erlang_d85+n_erlang_d15+1) = mu;
    end
    Rcycle(1+n_erlang_d85+n_erlang_d15,nstates-ncycles+cycle) = n_erlang_d15/y;
    R = [R; Rcycle];
end
R(end+1,end) = lambda;
R(end+1,:) = 0;

Q = (R-diag(R*ones(nstates,1)));
end