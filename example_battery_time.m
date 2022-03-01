function [nstates, npar, timer_aca, timer_reference, err] = example_battery_percentage(tol, ncycles, n_erlang_d85, n_erlang_d15, nt, np)


%% Model fixed parameters
% tol = 1e-6;
% ncycles = 2; 
lambda = 1e-5;

% approximate deterministic distribution with CDF delta(t-D85)
% as an Erlang with n_erlang_d85 states,
% and similarly for  the CDF delta(t-D15)
% n_erlang_d85 = 5;
% n_erlang_d15 = 3;

nstates = (ncycles-1)*(1+n_erlang_d85+n_erlang_d15+1)+1+1;

pi0 = zeros(nstates,1);
pi0(1) = 1;

% % y1 is DischargeTime and y2 is mu
% Q = @(y1,y2) infgen(nstates, ncycles, lambda, y1, y2, n_erlang_d85, n_erlang_d15);

%% Model variable parameters (time, DischargeTime and mu)
tf = 24*365*10;
p1low = 1;
p1up = 3;
p2low = 0.5;
p2up = 2.5;

% nt = 100; np = 100;

t = linspace(0, tf, nt);

DischargeTime = linspace(p1low, p1up, np);
mu = linspace(p2low, p2up, np);

n = [ length(t), length(DischargeTime), length(mu) ];
intervals = { t, DischargeTime, mu };

npar = prod(n);

% The infinitesimal generator as a function of the parameters \theta_1,
% ..., \theta_p.
Q = @(theta1, theta2) infgen(nstates, ncycles, lambda, theta1, theta2, n_erlang_d85, n_erlang_d15);

%% Percentage in D<15

kind = 'accumulated';
r = zeros(nstates,1);
for cycle = 1 : ncycles-1
    base_index = (cycle-1)*(1+n_erlang_d85+n_erlang_d15+1);
    for i = 0 : n_erlang_d15-1
        r(base_index+1+n_erlang_d85+i) = 1;
    end
end

% We create the function that evaluate a fiber of the tensor. 
Afiber = create_fiber_functions(Q, intervals, pi0, r, tol, kind);

% ACA
timer_aca = tic;
U = aca_nd(n, Afiber, tol);
timer_aca = toc(timer_aca);

% Reference solution
timer_reference = tic;
[RR, err] = create_reference_approximation(Afiber, intervals, U);
timer_reference = toc(timer_reference);

fprintf('Time, the full tensot has %d entries, the model has %d states\n', npar, nstates);

fprintf('Time for ACA: %f seconds\n', timer_aca);
fprintf('Time for explicit computation: %f seconds\n', timer_reference);
fprintf('Acceleration factor: %2.3fx\n', timer_reference / timer_aca);
fprintf('Relative error from the reference solution: %e\n', err);

save(strcat('case2_time_not_in_nr_RR_nstates_',int2str(nstates),'_nt_',int2str(nt),'_np_', int2str(np),'.mat'), "RR");

% plot Reliability at time tend
tend = floor(nt/5);
grid(:,:)=RR(tend,:,:)/(24*365*2);
f = figure;
imagesc(DischargeTime, mu, grid);
colorbar;
title('Percentage of time in 0<B$\le 15\%$ for $t_{\mathrm{end}}=2$ years','interpreter','latex');
xlabel('$L$','interpreter','latex'); 
ylabel('$\mu$','interpreter','latex'); 
saveas(f,strcat('case2_time_not_in_nr_nstates_',int2str(nstates),'_nt_',int2str(nt),'_np_', int2str(np),'.png'));

end

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
