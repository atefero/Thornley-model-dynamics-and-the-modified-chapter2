function [t,X] = simulateThornley(X0,t_end,params)
% This function simulates Thornley's model for the given timeframe over
% which to simulate 't_end', initial condition 'X0', and parameter values
% as specified in input struct 'params'
odeoptions = struct('AbsTol', 1e-6, 'RelTol', 1e-6);
[t,X] = ode15s( @(t,X) thornleyRHS(t,X,params), [0 t_end], X0, odeoptions);

end