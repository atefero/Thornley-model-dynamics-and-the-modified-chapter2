function params = applyThetaMod(params, theta_mod)
% This function takes the input set of parameters, and modifies the value
% of theta by the constant 'theta_mod'. All four parameters comprising
% theta are modified together, by portioning out the theta modifier into
% four parts.
modval = theta_mod^(1/4);
params.fN = params.fN * modval;
params.kC = params.kC * modval;
params.fC = params.fC / modval;
params.kN = params.kN / modval;

end