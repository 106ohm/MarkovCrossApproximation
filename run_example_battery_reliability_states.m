tol=1.e-6;
n_erlang_d85=10;
n_erlang_d15=5;
nt=10;
np=100;

for ncycles = 2 : 2 : 20 
    fprintf('RUN: ncycles=%d\n', ncycles);
    [nstates, npar, timer_aca, timer_reference, err] = example_battery_reliability(tol, ncycles, n_erlang_d85, n_erlang_d15, nt, np);               
end