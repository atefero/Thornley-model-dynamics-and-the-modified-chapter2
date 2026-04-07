function dX = reducedRHS(X,params)
% This function evaluates a reduced version of the right hand side of 
% Thornley's model, where the first two equations (which will be zero at a
% steady state point) are removed from the latter four equations. The
% species are still:
%   Msh - mass of dry matter in the plant shoot
%   Mrt - mass of dry matter in the plant root
%   Csh - concentration of substrate C in the plant shoot
%   Crt - concentration of substrate C in the plant root
%   Nsh - concentration of substrate N in the plant shoot
%   Nrt - concentration of substrate N in the plant root

% Growth/assimilation rates
Gsh = params.kG * X(3) * X(5);
Grt = params.kG * X(4) * X(6);
I = [ Gsh * X(1);
      Grt * X(2);
      params.kC / ( (1 + params.sigmaSS * X(1) / params.KM ) * (1 + params.sigmaPI * X(3) / params.JC ) );
      0;
      0;
      params.kN / ( (1 + params.sigmaSS * X(2) / params.KM ) * (1 + params.sigmaPI * X(6) / params.JN ) )   ];

% Litter rates - modified model adds litter terms to substrate C/N
Lsh = params.sigmaSS * params.klitt * X(1) / ( X(1) + params.KMlitt );
Lrt = params.sigmaSS * params.klitt * X(2) / ( X(2) + params.KMlitt );
L = [ -Lsh * X(1);
      -Lrt * X(2);
      -params.sigmaM * Lsh * X(3);
      -params.sigmaM * Lrt * X(4); 
      -params.sigmaM * Lsh * X(5); 
      -params.sigmaM * Lrt * X(6)             ];

% Consumption rates (due to growth)
C = [ 0;
      0;
      -params.fC * Gsh;
      -params.fC * Grt;
      -params.fN * Gsh;
      -params.fN * Grt    ];
      
% Transport rates
deltaC = X(3) - X(4);
deltaN = X(5) - X(6);
Afactor = 1 / X(1)^params.q + 1 / X(2)^params.q;
T = [ 0;
      0;
      -deltaC / ( X(1) * params.rhoC * Afactor );
       deltaC / ( X(2) * params.rhoC * Afactor );
      -deltaN / ( X(1) * params.rhoN * Afactor );
       deltaN / ( X(2) * params.rhoN * Afactor )   ];

% Combine the terms to get overall rate of change
dX = I + L + C + T;
         
end