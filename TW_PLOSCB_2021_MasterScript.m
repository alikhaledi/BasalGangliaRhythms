%% Master Script for Simulations/Figures in West et al. (2020) "State    %%
%%  Dependency of Cortical-Basal Ganglia Beta Rhythms"                   %%
% This is the main script that will reproduce the figures in the          %
% manuscript. This analysis takes forward a model of the cortico-basal    %
% ganglia circuit implemented in Dynamic Causal Modelling and previously  %
% described in van Wijk et al. (2018). This model was fit to data from    %
% West al. (2018) using an optimization based upon Approximate Bayesian   %
% Computation and described in West et al. (2021). Many of these analyses %
% are based upon scripts included in publicly available toolboxes and     %
% generously provided by their respective authors. Please see             %
% 'ABC_dependencies' for their respective licenses.                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Timothy West, Nuffield Department of Clinical Neurosciences, University %
% of Oxford; Wellcome Centre for Human Neuroscience, University College   %
% London. 2018-2021
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% TO DO:
% (1) Normalize spectra by state noise
restoredefaultpath
set(0,'defaultAxesFontSize',16)

clear; close all

% Initialise Analysis
%organise file paths, and add dependencies to path. Initializes
%configuration structure 'R' used throughout the scripts
R = project_AddPaths();
% Add configurations/settings to R
R = setupBasalGangliaModel(R);

% Figure (1): Model fit and comparison with data
simlength = 256; % Set simulation length
modID = 10; % This is the index of the winning model
BAA_plotModelFit(R,modID,simlength); %NOTE - You can dispense most of the wrappers - go straight to simulating data
plotDataComparison % plots the time series in figure 1A and C

%% Perform simulations with connectivity Sweep
BAA_sim_ConnectionSweep_v3(R)

% Figure (2): Sweep over connections and plot spectra
R = plotSweepSpectraBasic(R); % This doesnt just plot but also works out the vector Krange for spread of viable connections

% Computes Bursts from Connection Sweep Simulations
computeBurstWrapper_V3(R)
% Compute bursts from empirical data
computeBurstWrapper_V3_empiricalData(R)

% Model of Stimulation
BAA_sim_phaseLockedStim(R)
BAA_sim_phaseLockedStim_SICompEps(R)
peValidation
% Plot stim types
 BAA_plotStimTypeOutcomes(R)
 
% State matching analysis of stimulations
computeStimAnalysis_sweep(R,1); % the sweeps now do everything

%% Compare the connectivity matrices
fingerprintCompare_Sweep(R)% the sweeps now do everything
fingerprintCompare_Sweep_MSEALT(R); % with alternative cost functions
%% For Response Letter
computeOutOfBurstWrapper_V3(R); % compute spectral fingerprints without burst definition 
fingerprintCompareSimulatedEmpirical(R); %compare states


