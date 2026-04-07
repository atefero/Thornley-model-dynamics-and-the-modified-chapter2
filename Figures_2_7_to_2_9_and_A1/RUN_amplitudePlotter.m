function RUN_amplitudePlotter

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
params.sigmaM = 1;              % Selects original/modified model

% Specify the plotting colours
colors = [1.0, 0.0, 1.0;
          0.0, 0.0, 0.0;
          0.0, 0.6, 0.0;
          0.0, 0.0, 1.0;
          1.0, 0.0, 0.0];

% Specify the range of q values over which to plot, and how many values 
q_min = 0;
q_max = 8;
N_vals = 101;

% Specify the time range over which to pre-simulate
t_presim = 1e5;
% Specify the time over which to measure amplitudes after pre-simulation
t_measure = 1e3;

% Specify the initialisation values for the start of each loop
X0 = [2.417, 0.5246, 0.0246, 0.0090, 0.0084, 0.0142]';

% Specify the theta value for the first plot
plot_theta1 = 1;
% Specify the theta values to use for the second plot
plot_thetas2 = [0.1, 0.25, 0.5, 1.0, 1.5];

% Prepare ODE simulation options
odeoptions = struct('AbsTol', 1e-6, 'RelTol', 1e-6);

% Prepare q values to evaluate over, and storage for results
q_vec = linspace(q_min, q_max, N_vals);


%%% FIGURE ONE: PLOT THE BASE PARAMETER VALUES FIGURE

% Modify the the parameters to obtain the desired value of theta
params_here = applyThetaMod( params, plot_theta1 / (params.fN * params.kC / params.fC / params.kN) );

% Loop over all the q values
Mmin = zeros(1,N_vals);
Mmean = zeros(1,N_vals);
Mmax = zeros(1,N_vals);
for k = 1:N_vals
        
    % Set the q value to the current value in the loop
    params_here.q = q_vec(k);
    
    % Run the Thornley model for these parameters
    [T,Xpre] = ode15s(@(t,X) thornleyRHS(t,X,params_here), 0:t_presim, X0, odeoptions);
    
    % Now run it again for a shorter time period over which amplitudes are
    % measured, using a fine grid for better resolution
    [T,Xmeas] = ode15s(@(t,X) thornleyRHS(t,X,params_here), linspace(0,t_measure,t_measure*20), Xpre(end,:)', odeoptions);
    
    % Convert ODE output to total mass of dry matter
    M = Xmeas(:,1) + Xmeas(:,2);
    
    % Store minimum, mean and maximum
    Mmin(k) = min(M);
    Mmax(k) = max(M);
    Mmean(k) = mean(M);   
    
end

% Plot the results
figure; hold on;
plot(q_vec, Mmean, 'k--', 'LineWidth', 3);
plot(q_vec, Mmin, 'k', 'LineWidth', 3);
plot(q_vec, Mmax, 'k', 'LineWidth', 3);
% plot(q_vec, Mmax,'color', '#008000', 'LineWidth', 3);% , 'r' #31a354
% plot(q_vec, Mmean, 'k--', 'LineWidth', 3); % 'color', '#636363' ,
% plot(q_vec, Mmin, 'color', '#756bb1', 'LineWidth', 3);
%legend('Max M','Mean M','Min M','Location','best');
 xlabel('\bf{q}') % '\it\bf{q}'
 ylabel('\bf{M}')
 title('Modified');
 pubfig


 
%%% FIGURE TWO: PLOT THE DIFFERENT THETAS FIGURE

% Initialise the figure
figure; hold on;

% Prepare full storage 
N_thetas = length(plot_thetas2);
Mmin = zeros(N_thetas,N_vals);
Mmean = zeros(N_thetas,N_vals);
Mmax = zeros(N_thetas,N_vals);

% Loop over each theta value, and plot the results for its colour
parfor t = 1:N_thetas

    % Modify the the parameters to obtain the desired value of theta
    params_here = applyThetaMod( params, plot_thetas2(t) / (params.fN * params.kC / params.fC / params.kN) );
        
    % Loop over all the q values
    for k = 1:N_vals
    
        % Set the q value to the current value in the loop
        params_here.q = q_vec(k);
    
        % Run the Thornley model for these parameters
        [T,Xpre] = ode15s(@(t,X) thornleyRHS(t,X,params_here), 0:t_presim, X0, odeoptions);
    
        % Now run it again for a shorter time period over which amplitudes are
        % measured, using a fine grid for better resolution
        [T,Xmeas] = ode15s(@(t,X) thornleyRHS(t,X,params_here), linspace(0,t_measure,t_measure*20), Xpre(end,:)', odeoptions);
    
        % Convert ODE output to total mass of dry matter
        M = Xmeas(:,1) + Xmeas(:,2);
    
        % Store minimum, mean and maximum
        Mmin(t,k) = min(M);
        Mmax(t,k) = max(M);
        Mmean(t,k) = mean(M);   
    
    end
        
end

% Plot the results for this theta value, in its colour
for t = 1:N_thetas
    plot(q_vec, Mmean(t,:), '--', 'LineWidth', 3, 'Color', colors(t,:));
    plot(q_vec, Mmin(t,:), '-', 'LineWidth', 3, 'Color', colors(t,:));
    plot(q_vec, Mmax(t,:), '-', 'LineWidth', 3, 'Color', colors(t,:));
    legend( '' , '' , '\theta=0.1', '', '' , '\theta=0.25', '',  '', '\theta=0.5' , '' , '' , '\theta=1' , '',  '' , '\theta=1.5', 'FontSize', 14 ,'Location','best')
    xlabel('\bf{q}', 'fontsize' , 15)
    ylabel('\bf{M}', 'fontsize' , 15)
    title('Modified');
    hold on
    grid on
    pubfig
end