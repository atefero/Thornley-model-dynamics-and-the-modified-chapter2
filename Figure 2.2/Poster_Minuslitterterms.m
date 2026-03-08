function Poster_Minuslitterterms
% Generates dynamics and total-mass outputs for the modified Thornley model
% with litter terms applied to all six state variables.
%
% Associated with:
% "An in-depth study of the dynamics of Thornley's mathematical model in
% plant biology with a view to an improved model" (Rostami et al., 2025)
%
% Authors: Ati Rostami, Brodie A. J. Lawson, and Kevin Burrage

% 6 ODEs in the Modified Thornley Model:
    function u = fun_call(~,x)
        u = zeros(6,1);
        Q1 = x(3)/x(1) - x(4)/x(2);
        Q2 = pc/x(1)^q + pc/x(2)^q;

        R1 = x(6)/x(2) - x(5)/x(1);
        R2 = pn/x(2)^q + pn/x(1)^q;

        S = sigmaSS*klitt;

	%Shoot
        u(1) = kG * x(3) * x(5) / x(1) - S * x(1) / ( 1 + KMlitt / x(1) );
	
	%Root
        u(2) = kG * x(4) * x(6) / x(2) - S * x(2) /( 1 + KMlitt / x(2) );

	% Carbon Shoot Mass
        u(3) = kc * x(1) / ( ( 1 + sigmaSS * x(1) / KM ) * ( 1 + sigmaPI * x(3) / ( x(1) * Jc ) ) )...
            - fc * kG * x(3) * x(5) /x(1) - Q1/Q2 - S * x(3) / ( 1 + KMlitt / x(1) );

	% Carbon Root Mass
        u(4) = Q1/Q2 - fc * kG * x(4) * x(6) / x(2) - S * x(4) / ( 1 + KMlitt / x(2) );

	% Nitrogen Shoot Mass
        u(5) = R1/R2 - fn * kG * x(3) * x(5) / x(1) - S * x(5) / ( 1 + KMlitt / x(1) );

	% Nitrogen Root Mass
        u(6) = kn * x(2) / ( (1 + sigmaSS * x(2) / KM ) * (1 + sigmaPI * x(6)/( x(2) * Jn ) ) )...
            - fn * kG * x(4) * x(6) / x(2) - R1/R2 - S * x(6) / ( 1 + KMlitt / x(2) );
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
Mx = [0,0.9 ; 0, 2 ; 0,3.74 ; 0,4.32 ; 0.1,4.4178 ; 0,6];
%Mx = [0, 1; 0, 2;0,3.74; 0,4.32; 0.1,4.4179;0,6];
%Mx = [0, 0;0,2; 0,3.5; 0,4; 0,4.3;0,4.7];
%sigmaPI_vec =  [0, 0, 0, 0.1, 0.2, 0];                 % Product Inhibition   1(on)  0(off)
sigmaSS = 1;                        % Steady state growth
%q_vec = [2, 3.74, 4.29, 4.4179, 4.56, 6];                        % Change "transport resistance scaling parameter" to variable
%q_vec =6;
c = 0;                              % counter

for i = 1:size(Mx,1)      % Different Sigma PIs
    sigmaPI = Mx(i,1);       % Read out current sigmaPI from the vector
    q = Mx(i,2);             % Read out current q from the vector



    c = c+1;                    % Counter


    %Initial values defined in problem
    Msh0 = 0.8;
    Mrt0 = 0.9;
    Mshc0 = 0;
    Mrtc0 = 0.03;
    Mshn0 = 0.04;
    Mrtn0 = 0.001;

    x0 = [Msh0 , Mrt0 , Mshc0 , Mrtc0 , Mshn0 , Mrtn0 ];

    %% First solve from 0 to 20000
    tspan1 = 0:4000;
    x01 = x0;

    options = odeset('AbsTol',1e-6,'RelTol',1e-6);
    [t1,x1] = ode15s(@fun_call,tspan1,x01,options);

    mass_old = x1(:,1)+x1(:,2);


    %% Results for plotting

    M = (mass_old);               % For figure 3_Part A dry mass ; mass_new
    %D = diff(M(j))                   % Derivative of M

    % Figure 3 - Steady State
    figure (1)
    subplot(2,3,c);
     plot(t1, M,'Color','#2c7fb8','LineWidth', 2.5) 
    hold on
    %grid on
    %xlabel('\bf{Time} (Day)', 'fontsize' , 15)
    set (gca, 'fontsize', 15)
    %ylabel('\bf{Total Mass} (kg)', 'fontsize' , 15)
    title(['q = ',num2str(q),',   \sigma_{PI} = ',num2str(sigmaPI)],'FontSize',20);
    %title('Plant Structural Dry Mass')
    %legend({'\sigma_P_I = \it\bf{0}      \sigma_S_S = \it\bf{1}      q = \it\bf{[0, 2, 4.27, 6]}'},'FontSize', 20 , 'Location','best')
    legend( 'Original Model','Modified Model', 'FontSize', 14 ,'Location','best')

    pubfig

    %Plotting 6 State Variables:
    figure(3)
    subplot(2,3,c);
    plot(t1,x1(:,1),'b','linewidth',2);
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
    legend('Msh','Mrt','MshC','MrtC','MshN','MrtN', 'FontSize', 14 ,'Location','best');
    xlabel('\bf{Time}', 'fontsize' , 15)
    set (gca, 'fontsize', 15)
    ylabel('\bf{6 State Variables}', 'fontsize' , 15)
    %title('Modification of the Thornley model    6 ODEs')
    title(['q = ',num2str(q),',   \sigma_{PI} = ',num2str(sigmaPI)],'FontSize',20);

    drawnow
end
save test.mat

    figure(3)
    han=axes(figure(3),'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    xlabel(han,'Time (Day)','fontsize', 20); 
    ylabel(han,'6 State Variables','fontsize', 20);

    pubfig
end

