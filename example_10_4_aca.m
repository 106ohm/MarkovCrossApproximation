pi0 = [1.0;0;0;0];

tf = 10000.0;
plow = 0.5;
pup = 1.5;

t = linspace(0, tf, 100);
lambda1 = linspace(.5, 1.5, 100);
lambda2 = linspace(1e-4, 1e-3, 100);

Afiber = @(j,i) example_10_4_fiber(j, i, t, lambda1, lambda2, @Q, pi0);
n = [ length(t), length(lambda1), length(lambda2) ];
U = aca_nd(n, Afiber, 1e-6);

% pi5 = KolmogorovODE(Q([ 0.5 0.5 ]), pi0, t);
% pi15 = KolmogorovODE(Q([ 1.5 1.5 ]), pi0, t);

T = cpd2full(U);

% To check:
%
% pi10 = KolmogorovODE(Q_10_4([ lambda1(40), lambda2(60) ]), pi0, t);
% plot(t, 1-pi0(:,4), 'k*');
% hold on;
% plot(t, T(:,40,60));

function out = Q(y)
out = [-1-1e-3       1  1e-3    0
          0 -y(2)    0  y(2)
          0     0 -y(1) y(1)
          0     0    0    0];
end