f = @(p) p + sign(p - 0.35) * 0.2;
nparams = 4;
tol = 1e-2;

%% Model parameters
nreplicas = 4;

nstates = nreplicas+1; % labelled as n, n-1,..., 1 and 0 (failed system state)

pi0 = zeros(nstates,1);
pi0(1) = 1;

% theta is a vector of parameters, length(theta)=nparams, passed 
% entry by entry
density = 0.8;
R2 = sprand(nreplicas+1,nreplicas+1,density);
R2 = R2-diag(diag(R2));
Q = @(varargin)infgen(R2, f, nreplicas, varargin{:});

%% Model variable parameters (lambda and c)
tf = 10*nreplicas;
low = 0.25;
up = 0.5;

npoints = 50;
t = linspace(0, tf, npoints);
% lambda = repmat(linspace(low, up, npoints),nparams,1);
n = [ length(t), npoints*ones(1,nparams)];

intervals = {};
intervals{1} = t;
for p = 1 : nparams
    intervals{p+1} = linspace(low, up, npoints);
end

npar = prod(n);

%% Reliability at time t
kind = 'instantaneous';
r = [ones(nreplicas,1);0];

[Afiber, Aelem] = create_fiber_functions(Q, intervals, pi0, r, tol, kind);

% ACA
timer_aca = tic;
U = aca_nd(n, Afiber, tol);
timer_aca = toc(timer_aca);

%Reference solution
[RR, err] = create_reference_approximation(Afiber, intervals, U);

save('reference_50_p4.mat', "RR");