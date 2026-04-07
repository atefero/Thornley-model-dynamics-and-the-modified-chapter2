function J = thornleyJacobian(X,params)
% Author: Brodie A. J. Lawson
% This function evaluates the Jacobian matrix associated with the Thornley
% model, at the provided point X, for values of the parameters as specified
% in the input struct 'params'

% Initialise Jacobian (six equations, six variables)
J = zeros(6);
% Partial derivatives of RHS equation for mass of plant shoot
J(1,1) = params.kG * X(3) * X(5) - params.sigmaSS * params.klitt * X(1) * (X(1) + 2*params.KMlitt) / (X(1) + params.KMlitt)^2;
J(1,3) = params.kG * X(1) * X(5);
J(1,5) = params.kG * X(1) * X(3);
% Partial derivatives of RHS equation for mass of plant root
J(2,2) = params.kG * X(4) * X(6) - params.sigmaSS * params.klitt * X(2) * (X(2) + 2*params.KMlitt) / (X(2) + params.KMlitt)^2;
J(2,4) = params.kG * X(2) * X(6);
J(2,6) = params.kG * X(2) * X(4);
% Partial derivatives of RHS equation for concentration of substrate carbon in plant shoot
J(3,1) = -params.sigmaSS * params.kC / params.KM / (1 + params.sigmaSS * X(1)/params.KM)^2 / (1 + params.sigmaPI * X(3)/params.JC) + params.sigmaSS * (1 - params.sigmaM) * params.klitt * params.KMlitt * X(3) / (params.KMlitt + X(1))^2 + (X(3) - X(4)) * ( X(2)^(-params.q) + (-params.q+1) * X(1)^(-params.q)) / params.rhoC / ( X(1) * X(2)^(-params.q) + X(1)^(-params.q+1))^2;
J(3,2) = -(X(3) - X(4)) / params.rhoC / X(1) / ( X(1)^(-params.q) + X(2)^(-params.q) )^2 * params.q * X(2)^(-(params.q+1));
J(3,3) = -params.sigmaPI * params.kC / params.JC / (1 + params.sigmaSS * X(1)/params.KM) / (1 + params.sigmaPI * X(3)/params.JC)^2 + params.sigmaSS * (1 - params.sigmaM) * params.klitt * X(1) / (X(1) + params.KMlitt) - params.kG * X(5) * (params.fC + 2*X(3)) - 1/params.rhoC / X(1) / ( X(1)^(-params.q) + X(2)^(-params.q) );
J(3,4) = 1/params.rhoC / X(1) / ( X(1)^(-params.q) + X(2)^(-params.q) );
J(3,5) = -params.kG * X(3) * (X(3) + params.fC);
% Partial derivatives of RHS equation for concentration of substrate carbon in plant root
J(4,1) = (X(3) - X(4)) / params.rhoC / X(2) / ( X(1)^(-params.q) + X(2)^(-params.q) )^2 * params.q * X(1)^(-(params.q+1));
J(4,2) = params.sigmaSS * (1 - params.sigmaM) * params.klitt * params.KMlitt * X(4) / (params.KMlitt + X(2))^2 - (X(3) - X(4)) * ( X(1)^(-params.q) + (-params.q+1) * X(2)^(-params.q)) / params.rhoC / ( X(2) * X(1)^(-params.q) + X(2)^(-params.q+1))^2;
J(4,3) = 1/params.rhoC / X(2) / ( X(1)^(-params.q) + X(2)^(-params.q) );
J(4,4) = params.sigmaSS * (1 - params.sigmaM) * params.klitt * X(2) / (X(2) + params.KMlitt) - params.kG * X(6) * (params.fC + 2*X(4)) - 1/params.rhoC / X(2) / ( X(1)^(-params.q) + X(2)^(-params.q) );
J(4,6) = -params.kG * X(4) * (X(4) + params.fC);
% Partial derivatives of RHS equation for concentration of substrate nitrogen in plant shoot
J(5,1) = params.sigmaSS * (1 - params.sigmaM) * params.klitt * params.KMlitt * X(5) / (params.KMlitt + X(1))^2 + (X(5) - X(6)) * ( X(2)^(-params.q) + (-params.q+1) * X(1)^(-params.q)) / params.rhoN / ( X(1) * X(2)^(-params.q) + X(1)^(-params.q+1))^2;
J(5,2) = -(X(5) - X(6)) / params.rhoC / X(1) / ( X(1)^(-params.q) + X(2)^(-params.q) )^2 * params.q * X(2)^(-(params.q+1));
J(5,3) = -params.kG * X(5) * (X(5) + params.fN);
J(5,5) = params.sigmaSS * (1 - params.sigmaM) * params.klitt * X(1) / (X(1) + params.KMlitt) - params.kG * X(3) * (params.fN + 2*X(5)) - 1/params.rhoN / X(1) / ( X(1)^(-params.q) + X(2)^(-params.q) );
J(5,6) = 1/params.rhoN / X(1) / ( X(1)^(-params.q) + X(2)^(-params.q) );
% Partial derivatives of RHS equation for concentration of substrate nitrogen in plant root
J(6,1) = (X(5) - X(6)) / params.rhoC / X(2) / ( X(1)^(-params.q) + X(2)^(-params.q) )^2 * params.q * X(1)^(-(params.q+1));
J(6,2) = -params.sigmaSS * params.kN / params.KM / (1 + params.sigmaSS * X(2)/params.KM)^2 / (1 + params.sigmaPI * X(6)/params.JN) + params.sigmaSS * (1 - params.sigmaM) * params.klitt * params.KMlitt * X(6) / (params.KMlitt + X(2))^2 - (X(5) - X(6)) * ( X(1)^(-params.q) + (-params.q+1) * X(2)^(-params.q)) / params.rhoN / ( X(2) * X(1)^(-params.q) + X(2)^(-params.q+1))^2;
J(6,4) = -params.kG * X(6) * (X(6) + params.fN);
J(6,5) = 1/params.rhoN / X(2) / ( X(1)^(-params.q) + X(2)^(-params.q) );
J(6,6) = -params.sigmaPI * params.kN / params.JN / (1 + params.sigmaSS * X(2)/params.KM) / (1 + params.sigmaPI * X(6)/params.JN)^2 + params.sigmaSS * (1 - params.sigmaM) * params.klitt * X(2) / (X(2) + params.KMlitt) - params.kG * X(4) * (params.fN + 2*X(6)) - 1/params.rhoN / X(2) / ( X(1)^(-params.q) + X(2)^(-params.q) );

end