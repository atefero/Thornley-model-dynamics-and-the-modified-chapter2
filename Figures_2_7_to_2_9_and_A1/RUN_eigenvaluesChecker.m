function RUN_eigenvaluesChecker
%Figure 2.9
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
params.sigmaM = 0;              % Selects original/modified model

% Specify the maximum positive eigenvalue to plot
eig_max = 0.05;

% Specify the plotting colours
colors = [1.0, 0.2, 0.2;
          1.0, 0.2, 0.2;
          0.2, 0.2, 1.0;
          0.2, 0.2, 1.0];
      
% Specify the line styles (normal lines)
line_styles = {'-','--','-','--'};

% Specify the range of q values over which to plot, and how many values 
q_min = 0;
q_max = 6;
N_vals = 501;

% Set up the different switch variables to use across runs
sigmaPIs = [0,1,0,1];
sigmaMs = [0,0,1,1];

% Specify the theta values to use for plotting
plot_thetas = [0.25,1];
%plot_thetas = [0.75, 1.5];

% Prepare q values to evaluate over, and storage for results
q_vec = linspace(q_min, q_max, N_vals);
      
% Loop over different values of theta that were specified
for t = 1:length(plot_thetas)
    
    params_here = applyThetaMod( params, plot_thetas(t) / (params.fN * params.kC / params.fC / params.kN) );
    
    % Initialise Figure
    %figure; hold on;
        subplot(1,2,t);
        hold on
    
    % Add a marker for the zero line on the figure
    plot([q_min-1 q_max+1],[0 0],'k','LineWidth',3);
    
    fprintf('Working on theta value %g out of %g total\n\n',t,length(plot_thetas));
    
    for i = 1:length(sigmaPIs)
        
        fprintf(' ...working on switch setting %g out of %g total\n\n',i,length(sigmaPIs));
        
        % Use the loop to set product inhibition on or off
        params_here.sigmaPI = sigmaPIs(i);
        params_here.sigmaM = sigmaMs(i);
        
        % Prepare storage at beginning of loop
        valid = true(N_vals,1);
        Jeigs = zeros(N_vals,6);
        
        % Specify the initialisation values for the start of each loop
        X_star = [2.417, 0.5246, 0.0246, 0.0090, 0.0084, 0.0142]';
        
        for k = 1:N_vals
            
            % Set this to be the value of q
            params_here.q = q_vec(k);
            
            % Find the steady state values using the steady state finder function
            % (Use the previously found steady state as an initialisation)
            X_star = findSteadyState(params_here, true, X_star);
            
            % Check validity of steady state
            if any(X_star < 0) || any(~isreal(X_star))
                valid(k) = false;
            end
            
            % Determine the eigenvalues at steady state
            Jeigs(k,:) = eig(thornleyJacobian(X_star,params_here));
            
        end
        
        % Reduce the set of eigenvalues to only the eigenvalue with the
        % smallest absolute value real component
        [~,locs] = max( real(Jeigs), [], 2 );
        Jeigs = Jeigs(sub2ind([N_vals, 6],(1:N_vals)',locs));
        % Cap all very positive eigenvalues at a smaller number (for plot)
        % Includes a slight amount of jitter for display purposes
        jitter = eig_max / 20;
        Jeigs(Jeigs > eig_max) = eig_max + jitter*(i - 0.5 - length(sigmaPIs)/2);
        
        % Sort out the data by validity of steady state
        q_valid = nan(N_vals,1);
        q_valid(valid) = q_vec(valid);
        q_invalid = nan(N_vals,1);
        q_invalid(~valid) = q_vec(~valid);
        eigs_valid = nan(N_vals,1);
        eigs_valid(valid) = Jeigs(valid);
        eigs_invalid = nan(N_vals,1);
        eigs_invalid(~valid) = Jeigs(~valid);
        
        % Remove invalid data that is too close to other invalid data
        x_separation = 0.25;
        [q_invalid, eigs_invalid] = removeClose(q_invalid, eigs_invalid, x_separation, 1e-5);
        
        % Plot real components of Jacobian eigenvalues, coloured by validity
        plot(q_valid, real(eigs_valid), 'LineStyle', line_styles{i}, 'LineWidth', 2, 'Color', colors(i,:));
        plot(q_invalid, real(eigs_invalid), 'x', 'MarkerSize', 10, 'LineWidth', 3, 'Color', colors(i,:));
        
    end
    
    % Add the legend to the plot - plot on every even because two lines per
    % each run, and the first is the black line
    run_names = {'Original, \sigma_{PI} = 0', 'Original, \sigma_{PI} = 1', 'Modified, \sigma_{PI} = 0', 'Modified, \sigma_{PI} = 1'};
    for i = 1:length(run_names)
        legend_txts{2*i-1} = '';
        legend_txts{2*i} = run_names{i};
    end
    legend(legend_txts,'FontSize',24);
    
    % Label plot and set fontsize
    set(gca,'FontSize',24);
    xlabel('q','FontSize',24);
    ylabel('Max Re(\lambda)','Fontsize',24);
    % Adjust axis limits
    xlim([q_min q_max]);
    yl = ylim;
    ymin = yl(1);
    ymax = eig_max+(eig_max-ymin)*0.05;
    ylim([ymin ymax]);
    
    % Fix the ticks to reflect the capping
    yticks( linspace(ymin,eig_max,11) );
    yt_labels = yticklabels;
    yt_labels{end} = ['\geq',num2str(eig_max)];
    yticklabels(yt_labels);
    
    % Title shows value of theta
    title(['\theta = ',num2str(plot_thetas(t))],'FontSize',24);   
    pubfig
end

end

function [x,y] = removeClose(x,y,tx,ty)
% This function removes data points in (x,y) that are closer to the
% previous point in terms of x (distance tx) and y (distance ty). The 
% approach works sequentially along the input vectors x and y.

loc = 1;
while loc < length(x)
    
    % If the following point is too close, delete it
    if abs(x(loc+1) - x(loc)) < tx && abs(y(loc+1) - y(loc)) < ty
        x(loc+1) = [];
        y(loc+1) = [];
        
    % Otherwise, advance to the next location in the vector
    else
        loc = loc+1;
    end
    
end

end
    