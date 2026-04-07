function Surfforq(sigmaM)
% Author: Brodie A. J. Lawson
if sigmaM == 0
    filename = 'PARAMStest.xlsx';
else
    filename = 'modPARAMStest.xlsx';
end
sheet = 1;
qRange = 'A2:A133000';
sigmaRange ='B2:B133000';
Crange = 'C2:C133000';
q = xlsread(filename,sheet,qRange);
sigma_PI = xlsread(filename,sheet,sigmaRange);
C = xlsread(filename,sheet,Crange);

% Fill in with a colour per class for acustom color map
class_clr = [  17, 217, 65;           % Green
               7, 74, 242;           % Blue
               7, 172, 242;          % Light Blue
               242, 215, 10;          % yellow
               245, 133, 5;          % Orange
%168, 123, 50
             103, 108, 122;     % Grey
             ] / 256; % for example 5 color for 5 category

% Input the q levels used
q_vals = 0:0.025:6;
sigma_PI = 0:0.025:2;

% Create a grid of values using these
[Q, S] = meshgrid(q_vals, sigma_PI);

% Reshape the class data into a matrix of same size
C = reshape(C, length(q_vals), length(sigma_PI))';

% Define a shift to all colours that allows objects to be placed on top of
% the surface plot
Cshift = 50;

surf(Q,S,C-Cshift);
view(2);
shading flat;
% Tell MATLAB to use your colormap  (%Brodie Added The following to the code:)
colormap(class_clr);

% Label axes (and increase fontsize)
xlabel('\bf{q}','FontSize',30)
ylabel('\bf{\sigma_P_I}','FontSize',40)
title('Modified')

% Change position of axes
% Position:   [xstart, ystart, xlength, ylength]
set(gca,'Position',[0.05, 0.1, 0.5, 0.85]);
axis square;

% Tell MATLAB to add the colorbar to the plot (and give it a name)
clbar_obj = colorbar;
% Use name to give the colorbar new properties
clbar_obj.Ticks = [1, 2, 3, 4, 5, 6] - Cshift;
clbar_obj.TickLabels = {'Steady State', 'Damped Oscillations', 'Sustained Oscillations', 'Period Doubling/Chaos', 'Oscillation into Collapse', 'Death'};
clbar_obj.FontSize = 20;
clbar_obj.Position = [0.55, 0.05, 0.015, 0.65];
% Set axis limits for colorbar (so that 1, 2, 3, 4, 5, 6 are in middle)
caxis([0.5 6.5] - Cshift);
set (gca, 'fontsize', 24)
%pubfig

%%% ADD INSET IF VISUALISING THE MODIFIED MODEL

if sigmaM == 1
    
    % Define the position and dimensions of where to draw the inset box
%     inset_x = 0.325;
%     inset_y = 0.8;
%     inset_dx = 0.3;
%     inset_dy = 0.05;
    inset_x = 0.585;
    inset_y = 0.8;
    inset_dx = 0.3;
    inset_dy = 0.05;
    
    % Specify the region in (q,sigmaPI) space to plot
    inset_q = [2.4, 4.2];
    inset_sigmaPI = [0, 0.075];
    
    % Draw the box showing where the inset was taking from
    rectangle('Position',[inset_q(1), inset_sigmaPI(1), inset_q(2)-inset_q(1), inset_sigmaPI(2)-inset_sigmaPI(1)], 'LineWidth', 3, 'EdgeColor', [0 0 0]);
    
    % Store axis before creating new one
    orig_ax = gca;
        
    % Create the inset axes at the specified position
    inset_ax = axes('Position',[inset_x, inset_y, inset_dx, inset_dy]);
    
    % Subset the surf data to only those falling within the ranges
    %subrange = Q >= inset_q(1) & Q <= inset_q(2) & S >= inset_sigmaPI(1) & S <= inset_sigmaPI(2);
    %Qsub = Q(subrange);
    %Ssub = S(subrange);
    %Csub = C(subrange);
    
    % Replot the results inside the inset box
    surf(inset_ax, Q, S, C-Cshift);
    xlim(inset_ax, inset_q);
    ylim(inset_ax, inset_sigmaPI);
    shading(inset_ax, 'flat');
    view(inset_ax, 2);
    
    % Clean up the second axes
    box(inset_ax, 'on');
    set(inset_ax, 'FontSize', 20, 'LineWidth', 5);
    xticks( inset_q(1):0.6:inset_q(2) );
    %pubfig
end
    


%title('q vs. \sigma_P_I')
%legend({'1 = Steady State, 2 =  Damped Oscilations,  3 = Sustained Oscillations,  4 = Oscillations goes near zero,  5= Oscillations suddenly collapsed,  6 = Died,  7 = Just 2 Oscillations then suddenly collapsed'})