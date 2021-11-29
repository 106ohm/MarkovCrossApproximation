%% Model fixed parameters
tol = 1e-8;
ncycles = 7; 
lambda = 1e-4;

% approximate deterministic distribution with CDF delta(t-D85)
% as an Erlang with n_erlang_d85 states,
% and similarly for  the CDF delta(t-D15)
n_erlang_d85 = 5;
n_erlang_d15 = 3;

nstates = (ncycles-1)*(1+n_erlang_d85+n_erlang_d15+1)+1+1;

pi0 = zeros(nstates,1);
pi0(1) = 1;

% y1 is DischargeTime and y2 is mu
Q = @(y1,y2)infgen(nstates, ncycles, lambda, y1, y2, n_erlang_d85, n_erlang_d15);

%% Model variable parameters (time, DischargeTime and mu)
tf = 5.0;
p1low = 1;
p1up = 3;
p2low = 0.3;
p2up = 0.7;
mu = p2low;

t = linspace(0, tf, 400);


% pi is a matrix whose rows are pi(t_i) for t_i in [0,tf]
% pi = KolmogorovODE(Q(DischargeTime(1), mu(1)), pi0, t);

% DischargeTime = linspace(p1low, p1up, 400);
% mu = linspace(p2low, p2up, 400);
DischargeTime = chebpts(400, [p1low p1up]);
mu = chebpts(400, [p2low p2up]);

%U = ChebopMarkovOneParameter(@Q,pi0,tf,plow,pup);

%% Measure of interest
r = [ones(nstates-ncycles,1);zeros(ncycles,1)];% Reliability at time t

%% ACA related definitions
[rr,pp,kk] = rational(min(14, 5 + ceil(-log10(tol))), 1);
QQ = @(t1, t2) infgen(nstates, ncycles, lambda, t1, t2, n_erlang_d85, n_erlang_d15);
Afiber = @(j,i) example_battery_two_params_fiber(j,i,t,DischargeTime, mu, pi0, QQ, r, rr, pp, kk);

n = [ length(t), length(DischargeTime), length(mu) ];

U = aca_nd(n, Afiber, tol);

% Construct the reference solution
err = 0;
RR = zeros(length(t), length(DischargeTime), length(mu));
for i = 1 : length(DischargeTime)
    for j = 1 : length(mu)
        RR(:, i, j) = Afiber(1, [1, i, j]);
        
        % Construct the approximate sol
        as = zeros(length(t), 1);
        for k = 1 : length(t)
            as(k) = 0;
            for kk = 1 : size(U{1}, 2)
                as(k) = as(k) + U{1}(k, kk) * U{2}(i, kk) * U{3}(j, kk);
            end
        end
        
        err = hypot(err, norm( RR(:, i, j) - as, 'fro' ) );
    end
end

err / norm(RR(:))

%% Infinitesimal Generator Matrix
function Q = infgen(nstates, ncycles, lambda, y1, y2, n_erlang_d85, n_erlang_d15)
% the last ncycles states are (absorbing) failure states

% y1 is DischargeTime, y2 is mu
D85 = 0.85*y1;
D15 = 0.15*y1;
mu = y2;

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
        Rcycle(1+n_erlang_d85+i,base_index+1+n_erlang_d85+i+1) = n_erlang_d15/D15;
        Rcycle(1+n_erlang_d85+i,base_index+1+n_erlang_d85+n_erlang_d15+1) = mu;
    end
    Rcycle(1+n_erlang_d85+n_erlang_d15,nstates-ncycles+cycle) = n_erlang_d15/D15;
    R = [R; Rcycle];
end
R(end+1,end) = lambda;
R(end+1,:) = 0;

Q = (R-diag(R*ones(nstates,1)));
end
