%% Model fixed parameters
tol = 1e-6;
nreplicas = 2;

c2 = 0.9;
mu1 = 0.5;
mu2 = 0.5;

nstates = nreplicas+2; % labelled as n, n-1,..., 1 and 0 (failed system state)

pi0 = zeros(nstates,1);
pi0(1) = 1;

% theta1 is lambda (i.e., failure rate) and theta2 is c (i.e., coverage)
Q = @(theta1,theta2)infgen(nreplicas, theta1, theta2, c2, mu1, mu2);

%% Model variable parameters (lambda and c)
tf = 24*365*10;
p1low = 1.e-6;
p1up = 1.e-5;
p2low = 0.9;
p2up = 0.99;

npoints = 100;
t = linspace(0, tf, npoints);
lambda = linspace(p1low, p1up, npoints);
c = linspace(p2low, p2up, npoints);

n = [ length(t), length(lambda), length(c) ];
intervals = { t, lambda, c };

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

save('case1_modified_reliability_RR.mat', "RR");

% Use the GLM
tGLM = tic;
npoints_glm = 100;
X = zeros(npoints_glm, 3);
y = zeros(npoints_glm, 1);
f1 = @(x) x; f2 = @(x) x / min(abs(lambda)); f3 = @(x) x;
for i = 1 : npoints_glm
    X(i, :) = randi(npoints, 1, 3);
    y(i) = Aelem(X(i,:));
    X(i, :) = [ f1(t(X(i,1))), f2(lambda(X(i,2))), f3(c(X(i,3))) ];
end

mdl = fitglm(X, y, 'Distribution', 'Gamma'); tGLM = toc(tGLM);
% mdl = fitglm(X, y, 'poly233', 'Distribution', 'Gamma'); tGLM = toc(tGLM);
%plot(feval(mdl, X(:,1), X(:,2), X(:,3)), 'r-');
%hold on; plot(y, 'b-');

[xx,yy,zz] = ndgrid(1:length(t), 1:length(lambda), 1:length(c));
y = feval(mdl, f1(t(xx)), f2(lambda(yy)), f3(c(zz)));

errGLM = norm(RR(:) - y(:));
fprintf('errGLM = %e (time = %f)\n', errGLM / norm(RR(:)), tGLM);

norm(RR(:) - mean(RR(:))) / norm(RR(:))

% plot Reliability at time tend
tend = round(npoints/5);
grid = zeros(npoints, npoints);
grid(:,:)=RR(tend,:,:);
imagesc(lambda,c,grid);
colorbar;
title('Reliability at $t_{\mathrm{end}}=2$ years','interpreter','latex');
xlabel('$\lambda$','interpreter','latex'); 
ylabel('c'); 

%% Expected time under repair

kind = 'accumulated';
%r = zeros(nreplicas+1,1); r(1)=1;
r = ones(nreplicas+2,1);
r(1)=0; r(nreplicas+1)=0; r(nreplicas+2)=0;
 
Afiber = create_fiber_functions(Q, intervals, pi0, r, tol, kind);

% ACA
timer_aca = tic;
U = aca_nd(n, Afiber, tol);
timer_aca = toc(timer_aca);

% Reference solution
timer_reference = tic;
[RR, err] = create_reference_approximation(Afiber, intervals, U);
timer_reference = toc(timer_reference);

fprintf('Expected time not in nr, the full tensot has %d entries\n', npar);

fprintf('Time for ACA: %f seconds\n', timer_aca);
fprintf('Time for explicit computation: %f seconds\n', timer_reference);
fprintf('Acceleration factor: %2.3fx\n', timer_reference / timer_aca);
fprintf('Relative error from the reference solution: %e\n', err);

save('case1_modified_time_under_repair_RR.mat', "RR");

% plot 
tend = 20;
grid(:,:)=RR(tend,:,:);
imagesc(lambda,c,grid);
colorbar;
title('Expected time in $n_r-1,\dots,1$ with $t_{\mathrm{end}}=2$ years','interpreter','latex');
xlabel('$\lambda$','interpreter','latex'); 
ylabel('c'); 


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
