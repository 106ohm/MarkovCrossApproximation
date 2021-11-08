%% Model fixed parameters
ncycles = 7;
mu= 0.5; 
lambda = 1e-4;
D85 = 2;
D15 = 0.1;
tol = 1e-8;

% approximate deterministic distribution with CDF delta(t-D85)
% as an Erlang with n_erlang_d85 states,
% and similarly for  the CDF delta(t-D15)
n_erlang_d85 = 5;
n_erlang_d15 = 5;

nstates = (ncycles-1)*(1+n_erlang_d85+n_erlang_d15+1)+1+1;

pi0 = zeros(nstates,1);
pi0(1) = 1;

Q = @(y) infgen(nstates, ncycles, mu, lambda, D85, y, n_erlang_d85, n_erlang_d15);

%% Model variable parameters (time and D15)
tf = 5.0;
plow = 0.01;
pup = 0.5;

t = linspace(0, tf, 6000);

D15 = 0.1;
% pi is a matrix whose rows are pi(t_i) for t_i in [0,tf]
mpi = KolmogorovODE(Q(D15), pi0, t);

nd = 6000;
D15 = linspace(plow, pup, nd);

%U = ChebopMarkovOneParameter(@Q,pi0,tf,plow,pup);


%% Measure of interest
r = [ones(nstates-ncycles,1);zeros(ncycles,1)];% Reliability at time t

[rr,pp,kk] = rational(min(14, 5 + ceil(-log10(tol))), 1);

% Determine the permutation of Q that we should use
perm = amd(Q(D15(1)));
iperm = perm; for j = 1 : length(perm); iperm(perm(j)) = j; end

%% ACA related definitions
% Aelem = @(i,j) dot(r, KolmogorovPoint(Q(D15(j)), pi0, t(i)));
Aelem = @(i,j) dot(r, KolmogorovPoint2(Q(D15(j)), pi0, t(i), rr, pp, kk, perm, iperm));
Arowold  = @(i) arrayfun(@(j) Aelem(i, j), 1 : length(D15))';
Arow = @(i) rowChebWrapper(@(t, x) dot(r, KolmogorovPoint2(Q(x), pi0, t, rr, pp, kk, perm, iperm)), t(i), D15)';
Acol  = @(j) KolmogorovODE(Q(D15(j)), pi0, t) * r;

timer = tic;
[U, V] = aca_2d(length(t), length(D15), Aelem, Arow, Acol, tol);
toc(timer)

tic;
M = zeros(length(t), length(D15));
for j = 1 : length(D15)
    M(:,j) = Acol(j);
end
toc

norm(M - U*V', 'fro') / norm(M, 'fro')

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

Q = (R-spdiags(R*ones(nstates,1), 0, nstates, nstates));
end