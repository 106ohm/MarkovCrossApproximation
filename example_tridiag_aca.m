
m = 128;
pi0 = [1.0 ; zeros(m-1, 1) ];

tf = 100000.0;

plow = 1e-2;
pup = 1e-1;

t = linspace(0, tf, 100);

lambda1 = linspace(plow, pup, 100);
lambda2 = linspace(0.5, 1.5, 100);

Afiber = @(j,i) example_tridiag_fiber(j, i, t, lambda1, lambda2, m, @Q_tridiag, pi0);
n = [ length(t), length(lambda1), length(lambda2) ];
U = aca_nd(n, Afiber, 1e-3);

% pi5 = KolmogorovODE(Q([ 0.5 0.5 ]), pi0, t);
% pi15 = KolmogorovODE(Q([ 1.5 1.5 ]), pi0, t);

T = cpd2full(U);

% To check:
%
% pi10 = KolmogorovODE(Q_tridiag([ lambda1(40), lambda2(60) ], m), pi0, t);
% plot(t, 1-pi10(:,m), 'k*');
% hold on;
% plot(t, T(:,40,60));