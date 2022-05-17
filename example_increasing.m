%% variable number of parameters
nparams = 3;

%% Model parameters
tol = 1e-6;
nreplicas = 4;

nstates = nreplicas+1; % labelled as n, n-1,..., 1 and 0 (failed system state)

pi0 = zeros(nstates,1);
pi0(1) = 1;

% theta is a vector of parameters, length(theta)=nparams, passed 
% entry by entry
Q = @(varargin)infgen(nreplicas, varargin{:});

%% Model variable parameters (lambda and c)
tf = 10*nreplicas;
low = 0.25;
up = 0.5;

npoints = 2;
t = linspace(0, tf, npoints);
% lambda = repmat(linspace(low, up, npoints),nparams,1);
n = [ length(t), npoints*ones(1,nparams)];

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

% Reference solution
timer_reference = tic;
[RR, err] = create_reference_approximation(Afiber, intervals, U);
timer_reference = toc(timer_reference);

fprintf('Reliability, the full tensor has %d entries\n', npar);

fprintf('Time for ACA: %f seconds\n', timer_aca);
fprintf('Time for explicit computation: %f seconds\n', timer_reference);
fprintf('Acceleration factor: %2.3fx\n', timer_reference / timer_aca);
fprintf('Relative error from the reference solution: %e\n', err);

save(strcat('example_increasing_RR_',int2str(nparams),'.mat'), "RR");

norm(RR(:) - mean(RR(:))) / norm(RR(:))

% % plot Reliability at time tend
% tend = round(npoints/5);
% grid = zeros(npoints, npoints);
% grid(:,:)=RR(tend,:,:);
% imagesc(lambda,c,grid);
% colorbar;
% title('Reliability at $t_{\mathrm{end}}=2$ years','interpreter','latex');
% xlabel('$\lambda$','interpreter','latex'); 
% ylabel('c'); 

%% Infinitesimal Generator Matrix
function Q = infgen(nreplicas, varargin)
% birth-death process plus random deaths in group
br = 1;
dr = 1;
density = 0.3;

seed = 1;
rng(seed);

R1 = diag(dr*ones(nreplicas,1),1) + diag(br*ones(nreplicas,1),-1);

R2 = sprand(nreplicas+1,nreplicas+1,density);
R2 = R2-diag(diag(R2));
[i,j,~] = find(R2);
for k = 1 : size(varargin)
    R2(i(k),j(k)) = varargin{k};
end

R = R1 + R2;

Q = (R-diag(R*ones(nreplicas+1,1)));
end
