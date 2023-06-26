% This code was created by Loïc Marrec (loic.marrec@unibe.ch) and Thibault
% Bertrand (t.bertrand@imperial.ac.uk) - June 2023

clear all; close all; clc; 

GetSimulation = 'y';    % 'y' if you want to simulate data, 'n' otherwise
GetTheory = 'y';        % 'y' if you want to compute the probability, 'n' otherwise

% Global parameters
N0 = 1;             % Initial population size
K = 100;            % Carrying capacity
b = 1;              % Intrinsic birth rate 
GrowthModel = 'L';  % Growth Model: Blumberg ('B'), Gompertz ('G'), Logistic ('L'), Richards ('R')
if strcmp(GrowthModel,'B') || strcmp(GrowthModel,'R') 
    
    g = 1.5; % Additional parameter appearing in Blumberg and Richards models
    
else
    
    g = NaN;
    
end

tic;

% (1) Run Gillespie simulations 
if strcmp(GetSimulation, 'y')
    
    sim.T = 20;                                 % Last time point
    sim.Nt = 21;                                % Number of time points to save
    sim.t = linspace(0, sim.T, sim.Nt);         % List of time points
    sim.Nit = 1e3;                              % Number of stochastic realizations

    sim.N = GillespieStochasticGrowth(N0, K, b, g, sim.Nit, sim.t, GrowthModel);     % Simulation

    sim.PN = NaN(K, length(sim.t));             % Declare the table containing the simulated data
    
    for iN = 1 : K
        
        for iT = 1 : length(sim.t)
            
            sim.PN(iN, iT) = length(sim.N(iT, sim.N(iT, :) == iN))/sim.Nit;     % Compute the probability from the simulated data
            
        end
        
    end
    
    sim.PN = [zeros(1,length(sim.t)) ; sim.PN]; % Adding zero probability for N=0 case
          
    sim.Nmean = nanmean(sim.N(sim.N(:, length(sim.t)) ~= 0, :), 2);     % Compute the mean population size
    sim.Nstd = nanstd(sim.N(sim.N(:, length(sim.t)) ~= 0, :), 0, 2);    % Compute the standard deviation
    sim.Nci = 1.96.*sim.Nstd./sim.Nit;                                  % Compute the 95% confidence interval

    fname = ['PN_' GrowthModel '_sim.mat'];         % Save the data
    save(fname, 'N0', 'K', 'b', 'g', 'GrowthModel', 'sim');

    disp('--> Done with Gillespie simulations')
    toc
    
end

% (2) Compute the exact probability and the deterministic limit
if strcmp(GetTheory, 'y')
    
    th.T = 20;                        % Last time point
    th.Nt = 100;                      % Number of time points to save
    th.t = linspace(0, th.T, th.Nt);  % List of time points

    th.PN = stochasticGrowthPDF_serial(N0, K, b, g, th.t, GrowthModel);    % Compute the exact probability

    th.Nmean = sum(th.PN.*repmat((0 : K)', [1, length(th.t)]));   % Compute the mean population size

    fname = ['PN_' GrowthModel '_th.mat'];     % Save the data
    save(fname, 'N0', 'K', 'b', 'g', 'GrowthModel', 'th');
    disp('--> Done with computing the exact solution')
    toc
    
end

if strcmp(GetTheory, 'y') && strcmp(GetSimulation, 'y')

    fig = figure('Name', 'Results', 'NumberTitle', 'off');
    hold on  

        psim = errorbar(sim.t, sim.Nmean, sim.Nci, 'LineStyle', 'None', 'Marker', 'o', 'Color', 'r', 'MarkerFaceColor', 'r');
        pth = plot(th.t, th.Nmean, 'LineStyle', '-', 'Linewidth', 1.5, 'Marker', 'None', 'Color', 'b');

    hold off
    hXLabel = xlabel('Time t', 'Color', 'k');
    hYLabel = ylabel('Population size N', 'Color', 'k');
    hLegend = legend([psim pth], {'Simulation', 'Exact'});
    set( gca                       , ...
        'FontName'   , 'Arial'   , 'FontSize'   , 14);
    set([hXLabel, hYLabel], ...
        'FontName'   , 'Arial'   , 'FontSize'   , 14);
    set(hLegend, ...
        'FontName'   , 'Arial'   , 'FontSize'   , 12, 'Location', 'East');
    set(gca, ...
      'Box'         , 'off'     , ...
      'TickDir'     , 'out'     , ...
      'TickLength'  , [.02 .02] , ...
      'XMinorTick'  , 'on'      , ...
      'YMinorTick'  , 'on'      , ...
      'YGrid'       , 'off'      , ...
      'LineWidth'   , 1         );
    ax = gca;
    ax.XAxis.Color = 'k';
    ax.YAxis.Color = 'k';
    axis tight
    ylim([0 K])
    
end