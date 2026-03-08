function ZoomInq2
% Original (Thornley) Model: The Transport-Resistance Model of Allocation (Thornley 1997)
% Figure 2.4
% Generates total mass trajectories and six-state-variable dynamics
% for the original Thornley model under different values of q and sigma_PI.
%
% Authors: Ati Rostami, Brodie A. J. Lawson, and Kevin Burrage


% 6 ODEs in Thornley Paper:
    function u = fun_call(~,x)
        
        %x = [max([zeros(1,6);x'])]';
        
        u = zeros(6,1);
        Q1 = x(3)/x(1) - x(4)/x(2);
        Q2 = pc/x(1)^q + pc/x(2)^q;
        
        R1 = x(6)/x(2) - x(5)/x(1);
        R2 = pn/x(2)^q + pn/x(1)^q;
        
        S = sigmaSS*klitt;
        
        u(1) = kG * x(3) * x(5) / x(1) - S * x(1) / ( 1 + KMlitt / x(1) );
        
        u(2) = kG * x(4) * x(6) / x(2) - S * x(2) /( 1 + KMlitt / x(2) );
        
        u(3) = kc * x(1) / ( ( 1 + sigmaSS * x(1) / KM ) * ( 1 + sigmaPI * x(3) / ( x(1) * Jc ) ) )...
            - fc * kG * x(3) * x(5) /x(1) - Q1/Q2;
        
        u(4) = Q1/Q2 - fc * kG * x(4) * x(6) / x(2);
        
        u(5) = R1/R2 - fn * kG * x(3) * x(5) / x(1);
        
        u(6) = kn * x(2) / ( (1 + sigmaSS * x(2) / KM ) * (1 + sigmaPI * x(6)/( x(2) * Jn ) ) )...
            - fn * kG * x(4) * x(6) / x(2) - R1/R2;
    end


%Define Constants
kG = 200;       % Growth Rate
klitt = 0.05;   % Litter Rate
KMlitt = 0.5;   % Litter Parameter
Jc = 0.1;       % Inhibition of C assimilation
Jn = 0.01;      % Inhibition of N uptake
fc = 0.5;       % Fractions of C in structural   %%fc = 0.5;
fn = 0.025;     %% fN 0.025 is better  % Fractions of N in structural %%fn = 0.025;
pc = 1;         % Transport resistance coefficients
pn = 1;         % Transport resistance coefficients
kc = 0.1;       % C assimilation parameter    %%kc = 0.1;
kn = 0.02;      % N uptake parameter          %%kn = 0.02;
KM = 1;         % Parameter giving asymptotic values of photosynthesis and N uptake



% Switches and q
tspan1 = 0:80000; %% time vector
sigmaPI_vec =  0;                 % Product Inhibition   1(on)  0(off)
sigmaSS = 1;                        % Steady state growth
q_vec = 2;                        % Change "transport resistance scaling parameter" to variable
c = 0;                              % counter
num_iterations = length(q_vec);
M_values = zeros(length(tspan1), num_iterations);

for i = 1:length(sigmaPI_vec)       % Different Sigma PIs
    sigmaPI = sigmaPI_vec(i);       % Read out current sigmaPI from the vector
    
    for j = 1:length(q_vec)         % different qs
        q = q_vec(j);               % Read out current q from the vector
        
        
        c = c+1;                    % Counter
        sigmaPI_list(c) = sigmaPI;    % Long list of all the different sigmaPIs for every single individual run
        q_list(c) = q;              % Long list of all the different qs for every single individual run
        
        
         %Initial values defined in problem
    Msh0 = 0.8;
    Mrt0 = 0.9;
    Mshc0 = 0;
    Mrtc0 = 0.03;
    Mshn0 = 0.04;
    Mrtn0 = 0.001;
    x0 = [Msh0, Mrt0, Mshc0, Mrtc0, Mshn0, Mrtn0];
    % Uncomment to run from apparent steady state
    % load('steady_State.mat','Msh_star','Mrt_star','MshC_star','MrtC_star','MshN_star','MrtN_star');
    % x0 = [Msh_star , Mrt_star , MshC_star , MrtC_star , MshN_star , MrtN_star ];
    
    %% First solve from 0 to 20000
    %tspan1 = 0:70000;
    x01 = x0;
    
    options = odeset('AbsTol',1e-10,'RelTol',1e-10);
    [t1,x1] = ode15s(@fun_call,tspan1,x01,options);
    
    mass_old = x1(:,1)+x1(:,2);
    
    
    %% Results for plotting
    
    M = (mass_old);               % For figure 3_Part A dry mass ; mass_new

        % Figure 3 - Steady State
        figure (1)
        plot(t1, M,'LineWidth', 2.5)
        %         [pks,locs, widths,proms] = findpeaks(M,t1)
        %         findpeaks(M,t1)
        %         text(locs+.02,pks,num2str((1:numel(pks))'))
        hold on
        grid on
        xlabel('\bf{Time}')
        ylabel('\bf{M}')
        title('Plant Structural Dry Mass')
        legend({'\sigma_P_I = \it\bf{0.4}      \sigma_S_S = \it\bf{1}      q = \it\bf{4.27}'}, 'FontSize', 20 , 'Location','best')
        %legend({'\sigma_P_I = 0  \sigma_S_S = 1  q = 3','\sigma_P_I = 0  \sigma_S_S = 1  q = 1','\sigma_P_I = 1  \sigma_S_S = 1  q = 2','\sigma_P_I = 1  \sigma_S_S = 1  q = 3','\sigma_P_I = 0  \sigma_S_S = 1  q = 4','\sigma_P_I = 0  \sigma_S_S = 1  q = 5'} , 'Location','best');
        
        
        %Plotting 6 State Variables:
        figure(3)
        plot(t1,x1(:,1),'b','linewidth',2.5);
        hold on
        plot(t1,x1(:,2),'k','linewidth',2);
        hold on
        plot(t1,x1(:,3),'g','linewidth',2);
        hold on
        plot(t1,x1(:,4),'r','linewidth',2);
        hold on
        plot(t1,x1(:,5),'c','linewidth',2);
        hold on
        plot(t1,x1(:,6),'m','linewidth',2);
        legend('Msh','','','','','', 'FontSize', 40 ,'Location','best');
        xlabel('\bf{Time}', 'fontsize' , 40)
        set (gca, 'fontsize', 40)
        ylabel('\bf{Msh}', 'fontsize' , 40)
        title(['q = ',num2str(q),',   \sigma_{PI} = ',num2str(sigmaPI)],'FontSize',30);
%           xlim([10000,70000]);
%           ylim([2.435,2.445]);
          xlim([14985,64221]);
          ylim([2.436,2.444]);
        %% 2-D scatter
        % Plot the grid of all the different combination
        %         figure (2)
        %         scatter(q_list,sigmaPI_list,50, [1] , 'filled');
        %         hold on
        %         grid on
        %         xlabel('q')
        %         ylabel('\sigma_P_I ')
        %         title('\sigma_P_I  vs. q')
        %         legend({'Plant is alive', 'Low peak Oscillation','High peak Oscillation','Died out'} , 'Location','best');
        %
        pubfig
        drawnow
    end
    save ZoomInq2.mat
end
end

