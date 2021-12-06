for tol = linspace(1.e-8,1.e-6,3)
    for ncycles = 2 : 5
        for n_erlang_d85 = 4 : 4 : 8
            for n_erlang_d15 = 2 : 2 : 4
                for nt = 10 : 10 : 100
                    for np = 2 : 2 : 10
                        fprintf('RUN: tol=%e, ncycles=%d, ne85=%d, ne15=%d, nt=%d, np=%d\n', tol, ncycles, n_erlang_d85, n_erlang_d15, nt, np);
                        [nstates, npar, timer_aca, timer_reference, err] = example_battery_reliability(tol, ncycles, n_erlang_d85, n_erlang_d15, nt, np);
                    end
                end
            end
        end
    end
end