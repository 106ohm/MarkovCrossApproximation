%% Model fixed parameters
tol = 1e-6;
nreplicas = 2;
mu = 0.5;

nstates = nreplicas+1; % labelled as n, n-1,..., 1 and 0 (failed system state)

pi0 = zeros(nstates,1);
pi0(1) = 1;

% theta1 is lambda (i.e., failure rate) and theta2 is c (i.e., coverage)
Q = @(theta1,theta2)infgen(nreplicas, theta1, theta2, mu);

%% Model variable parameters (lambda and c)
tf = 24*365*10;
p1low = 1.e-6;
p1up = 1.e-5;
p2low = 0.9;
p2up = 0.99;

npoints = 20;
t = linspace(0, tf, npoints);
lambda = linspace(p1low, p1up, npoints);
c = linspace(p2low, p2up, npoints);

n = [ length(t), length(lambda), length(c) ];
intervals = { t, lambda, c };

npar = prod(n);

%% Reliability at time t
kind = 'instantaneous';
r = [ones(nreplicas,1);0];

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

save('case1_reliability_RR.mat', "RR");

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

keyboard

% plot Reliability at time tend
tend = round(npoints/5);
grid = zeros(npoints, npoints);
grid(:,:)=RR(tend,:,:);
imagesc(lambda,c,grid);
colorbar;
title('Reliability at $t_{\mathrm{end}}=2$ years','interpreter','latex');
xlabel('$\lambda$','interpreter','latex'); 
ylabel('c'); 

%% Energy consumption

kind = 'accumulated';
%r = zeros(nreplicas+1,1); r(1)=1;
r = ones(nreplicas+1,1);
r(1)=0; r(nreplicas+1)=0;
 
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

save('case1_time__not_nr_RR.mat', "RR");

% plot 
tend = 20;
grid(:,:)=RR(tend,:,:);
imagesc(lambda,c,grid);
colorbar;
title('Expected time in $n_r-1,\dots,1$ with $t_{\mathrm{end}}=2$ years','interpreter','latex');
xlabel('$\lambda$','interpreter','latex'); 
ylabel('c'); 

% %% Energy consumption
% 
% kind = 'accumulated';
% r = linspace(nreplicas,0,nreplicas+1)';
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
% fprintf('Energy consumption, the full tensot has %d entries\n', npar);
% 
% fprintf('Time for ACA: %f seconds\n', timer_aca);
% fprintf('Time for explicit computation: %f seconds\n', timer_reference);
% fprintf('Acceleration factor: %2.3fx\n', timer_reference / timer_aca);
% fprintf('Relative error from the reference solution: %e\n', err);
% 
% save('case1_energy_RR.mat', "RR");
% 
% % plot Energy consumed at time tend
% tend = 20;
% Power = 10; %Watts
% grid(:,:)=Power*RR(tend,:,:)/1000;
% imagesc(lambda,c,grid);
% colorbar;
% title('Expected energy (kWh) consumed at $t_{\mathrm{end}}=2$ years','interpreter','latex');
% xlabel('$\lambda$','interpreter','latex'); 
% ylabel('c'); 

%% Infinitesimal Generator Matrix
function Q = infgen(nreplicas, y1, y2, mu)
% the working states are "nreplicas", "nreplicas-1", ..., 1,
% the last state is (absorbing) system failure state, here there are no
% working replicas

R = sparse(nreplicas+1,nreplicas+1);
for i = nreplicas : -1 : 2
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
