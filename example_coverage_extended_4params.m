%% Model fixed parameters
tol = 1e-6;
nreplicas = 2;
mu1 = 0.5;

nstates = nreplicas+2; % labelled as n, n-1,..., 1, standby and 0 (failed system state)

pi0 = zeros(nstates,1);
pi0(1) = 1;

% theta1 is lambda (i.e., failure rate), theta2 is c1 (i.e., failure coverage)
% theta3 is c2 (i.e., recovery coverage) and theta4 is mu2 (i.e., recovery rate)
Q = @(theta1,theta2,theta3,theta4)infgen(nreplicas, theta1, theta2, theta3, mu1, theta4);

%% Model variable parameters (lambda, c1 and c2)
tf = 24*365*10;
p1low = 1.e-6;
p1up = 1.e-5;
p2low = 0.9;
p2up = 0.99;
p3low = 0.9;
p3up = 0.99;
p4low = 0.25;
p4up = 0.75;

npoints = 16;
t = linspace(0, tf, npoints);
lambda = linspace(p1low, p1up, npoints);
c1 = linspace(p2low, p2up, npoints);
c2 = linspace(p3low, p3up, npoints);
mu2 = linspace(p4low, p4up, npoints);

n = [ length(t), length(lambda), length(c1), length(c2), length(mu2) ];
intervals = { t, lambda, c1, c2, mu2 };

npar = prod(n);

%% Reliability at time t
kind = 'instantaneous';
r = [ones(nreplicas+1,1);0];

[Afiber, Aelem] = create_fiber_functions(Q, intervals, pi0, r, tol, kind);

% ACA
timer_aca = tic;
U = aca_nd(n, Afiber, tol);
timer_aca = toc(timer_aca);

% Reference solution
timer_reference = tic;
[RR, err] = create_reference_approximation(Afiber, intervals, U);
timer_reference = toc(timer_reference);

fprintf('Reliability, the full tensor has %d entries\n', npar);

fprintf('Time for ACA: %f seconds\n', timer_aca);
fprintf('Time for explicit computation: %f seconds\n', timer_reference);
fprintf('Acceleration factor: %2.3fx\n', timer_reference / timer_aca);
fprintf('Relative error from the reference solution: %e\n', err);

save('case1_extended_reliability_RR.mat', "RR");

% %% Expected time under repair
% 
% kind = 'accumulated';
% %r = zeros(nreplicas+1,1); r(1)=1;
% r = ones(nreplicas+1,1);
% r(1)=0; r(nreplicas+1)=0;
%  
% Afiber = create_fiber_functions(Q, intervals, pi0, r, tol, kind);
% 
% % ACA
% timer_aca = tic;
% U = aca_nd(n, Afiber, tol);
% timer_aca = toc(timer_aca);
% 
% % Reference solution
% timer_reference = tic;
% [RR, err] = create_reference_approximation(Afiber, intervals, U);
% timer_reference = toc(timer_reference);
% 
% fprintf('Expected time under repair, the full tensot has %d entries\n', npar);
% 
% fprintf('Time for ACA: %f seconds\n', timer_aca);
% fprintf('Time for explicit computation: %f seconds\n', timer_reference);
% fprintf('Acceleration factor: %2.3fx\n', timer_reference / timer_aca);
% fprintf('Relative error from the reference solution: %e\n', err);
% 
% save('case1_extended_time_not_nr_RR.mat', "RR");

%% Infinitesimal Generator Matrix
function Q = infgen(nreplicas, lambda, c1, c2, mu1, mu2)
% the working states are "nreplicas", "nreplicas-1", ..., 1,
% the last minus one is a standby state (where the system has
% the last change to recover),
% and the last state (absorbing) is system failure state

R = sparse(nreplicas+2,nreplicas+2);
for i = nreplicas : -1 : 2
    % i counts the number of working replicas
    R(nreplicas-i+1, nreplicas-i+2) = i*lambda*c1;
    R(nreplicas-i+1,nreplicas+1) = i*lambda*(1-c1);
    if i<nreplicas
        R(nreplicas-i+1,nreplicas-i) = mu1;
    end
end
R(nreplicas, nreplicas+1) = lambda;
R(nreplicas, nreplicas-1) = mu1;

R(nreplicas+1,1) = c2*mu2;% rejuvenation
R(nreplicas+1,nreplicas+2) = (1-c2)*mu2;

Q = (R-diag(R*ones(nreplicas+2,1)));
end
