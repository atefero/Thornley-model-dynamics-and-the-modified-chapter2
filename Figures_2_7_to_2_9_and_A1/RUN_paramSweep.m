function RUN_paramSweep(sigmaM)
% This function takes as input sigmaM (0 for original model and 1 for
% modified model):   In Command Window:   RUN_paramSweep(0)  or RUN_paramSweep(1)
% Figure 2.7
% Authors: Ati Rostami, Brodie A. J. Lawson, and Kevin Burrage
% Define the parameter values (all but q, which will be varied)
params.kG = 200;                 % Growth rate
params.fC = 0.5;                 % Fraction of C in dry matter
params.fN = 0.025;               % Fraction of N in dry matter
params.KM = 1;                   % Saturation constant for C/N assimilation
params.kC = 0.1;                 % Assimilation rate for C
params.kN = 0.02;                % Uptake rate for N
params.JC = 0.1;                 % Inhibition rate for C (product inhib.)
params.JN = 0.01;                % Inhibition rate for N (product inhib.)
params.klitt = 0.05;             % Litter rate
params.KMlitt = 0.5;             % Saturation constant for litter
params.rhoC = 1;                 % Transport resistance
params.rhoN = 1;                 % Transport resistance

% Switches treated as parameters here
params.sigmaSS = 1;
params.sigmaPI = 0;
params.sigmaM = sigmaM;              % Selects original/modified model           

% Specify the range for the values to loop over
q_min = 0;
q_max = 6;
sigmaPI_min = 0;
sigmaPI_max = 2;
dp = 0.025;

% Specify the cutoffs used for automatic determination of phenomena
M_death = 1e-3;    % Total mass that must be exceeded for plant to be considered alive
Mdot_osc = 1e-3;   % Derivative of total mass that must be exceeded for plant to be considered oscillating
Pmin = 0.05;       % Minimum prominence for peaks to be registered
dPmax = 0.01;      % Difference in peak prominence tolerated for "equivalent" peaks

% Specify the time range over which to pre-simulate
t_presim = 1e5;
% Specify the time over which to measure amplitudes after pre-simulation
t_measure = 5e3;

% Specify the initialisation values for each simulation
X0 = [0.8, 0.9, 0, 0.03/0.9, 0.04/0.8, 0.001/0.9]';

% Prepare ODE simulation options
odeoptions = struct('AbsTol', 1e-6, 'RelTol', 1e-6);

% Prepare the grid of values to loop over
[P,Q] = meshgrid(sigmaPI_min:dp:sigmaPI_max, q_min:dp:q_max);
% Store these as vectors
Pv = P(:);
Qv = Q(:);

% Check number of sims that will be run
N_sims = length(Pv);

% Prepare result storage
Rv = zeros(N_sims,1);
oscillation_before = zeros(N_sims,1);
oscillation_after = zeros(N_sims,1);

% Loop over all the values of the parameters
for k = 1:N_sims
        
    % Set the parameters here
    params.sigmaPI = Pv(k);
    params.q = Qv(k);
    
    % Run the model for a very long time
    [Tpre,Xpre] = ode15s(@(t,X) thornleyRHS(t,X,params), 0:t_presim, X0, odeoptions);
    % Run it for a shorter, "measurement" period of time
    [Tmeas,Xmeas] = ode15s(@(t,X) thornleyRHS(t,X,params), linspace(0,t_measure,t_measure*20), Xpre(end,:)', odeoptions);
    
    % Calculate total mass versions of the model output
    Mpre = Xpre(:,1) + Xpre(:,2);
    Mmeas = Xmeas(:,1) + Xmeas(:,2);
    
    % Check for peaks in the initial stages, used for auto-classification
    pks = findpeaks(Mpre(Tpre < 0.1*t_presim),'MinPeakProminence',Pmin);
        
    % Check if the plant showed any life during the final window
    if any(Mmeas > M_death)
        
        % Check if there is movement in the late stages
        if any(abs(diff(Mmeas)./diff(Tmeas)) > Mdot_osc, 'all')
            
            % Check for peaks in the final period to examine the nature of
            % the oscillations - record their prominences
            [~,~,~,P] = findpeaks(Mmeas, 'minPeakProminence',Pmin);
            
            % Check if there are sufficient peaks to evaluate here
            if length(P) >= 5
                
                % Remove the first and last prominence as these will be
                % affected by cutoff at the boundaries
                P(1) = [];
                P(end) = [];
            
                % If the oscillations are sufficiently consistent, record
                % this as sustained oscillations
                if all( abs(P - mean(P)) < dPmax )          
                    Rv(k) = 3;         % SUSTAINED OSCILLATIONS
                else
                    Rv(k) = 4;         % CHAOTIC OSCILLATIONS
                end
                
            else
                Rv(k) = NaN;           % DON'T KNOW BUT AVOID ERROR
            end
            
        else
            
            % Check if there are any significant oscillations
            if length(pks) > 1
                Rv(k) = 2;     % DAMPED OSCILLATIONS
            else
                Rv(k) = 1;     % STEADY STATE
            end
        end
        
    else
        
        % Check if there is more than one peak, with growing amplitude
        if length(pks) > 1 && ~isequal(pks, sort(pks,'descend'))
            Rv(k) = 5;        % GROWING OSCILLATIONS INTO DEATH
        else
            Rv(k) = 6;
        end
        
    end
    
    % Output progress to user
    if mod(k,1000) == 0
        fprintf('Completed simulation %g out of %g\n', k, N_sims);
    end
    
end

% Create a table for file saving
Tmeas = table(Qv,Pv,Rv,ones(N_sims,1),oscillation_before,oscillation_after,'VariableNames',{'q','sigma_PI','C','sigma_SS','oscBefore','oscAfter'});

% Save the results in the Excel file
if sigmaM == 0
    writetable(Tmeas,'PARAMStest.xlsx');
else
    writetable(Tmeas,'modPARAMStest.xlsx');
end