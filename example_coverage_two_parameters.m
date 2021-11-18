%% Model fixed parameters
tol = 1e-6;
nreplicas = 5; 
mu = 0.5;

nstates = nreplicas+1; % labelled as n, n-1,..., 1 and 0 (failed system state)

pi0 = zeros(nstates,1);
pi0(1) = 1;

% y1 is lambda (i.e., failure rate) and y2 is c (i.e., coverage)
Q = @(y1,y2)infgen(nreplicas, y1, y2, mu);

%% Model variable parameters (lambda and c)
tf = 5.0;
p1low = 0.001;
p1up = 0.01;
% lambda = p1low;
p2low = 0.9;
p2up = 0.99;
% c = p2low;

t = linspace(0, tf, 100);

% pi is a matrix whose rows are pi(t_i) for t_i in [0,tf]
% pi = KolmogorovODE(Q(lambda, c), pi0, t);

lambda = linspace(p1low, p1up, 100);
c = linspace(p2low, p2up, 100);

%U = ChebopMarkovOneParameter(@Q,pi0,tf,plow,pup);

%% Measure of interest
r = [ones(nreplicas,1);0];% Reliability at time t

%% ACA related definitions
[rr,pp,kk] = rational(min(14, 5 + ceil(-log10(tol))), 1);
QQ = @(t1, t2) infgen(nreplicas, t1, t2, mu);
Afiber = @(j,i) example_coverage_two_params_fiber(j,i,t,lambda, c, pi0, QQ, r, rr, pp, kk);

n = [ length(t), length(lambda), length(c)];

U = aca_nd(n, Afiber, tol);

% Construct the reference solution
err = 0;
RR = zeros(length(t), length(lambda), length(c));
for i = 1 : length(lambda)
    for j = 1 : length(c)
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
function Q = infgen(nreplicas, y1, y2, mu)
% the working states are "nreplicas", "nreplicas-1", ..., 1,
% the last state is (absorbing) system failure state, here there are no
% working replicas

R = sparse(nreplicas+1,nreplicas+1);
for i = nreplicas : 2
    % i counts the number of working replicas
    R(nreplicas-i+1, nreplicas-i+2) = i*y1*y2;
    R(nreplicas-i+1,nreplicas+1) = i*y1*(1-y2);
    if i<nreplicas
        R(nreplicas-i+1,nreplicas-i) = mu;
    end
end
R(nreplicas, nreplicas+1) = y1;
R(nreplicas, nreplicas-1) = mu;

Q = (R-diag(R*ones(nreplicas+1,1)));
end
