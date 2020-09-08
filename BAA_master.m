%% Master Script for Simulations/Figures in West et al. (2020) "State    %%  
%%  Dependency of Cortical-Basal Ganglia Beta Rhythms"                   %%
% This is the main script that will reproduce the figures in the          %
% manuscript. This analysis takes forward a model of the cortico-basal    %
% ganglia circuit implemented in Dynamic Causal Modelling and previously  %
% described in van Wijk et al. (2018). This model was fit to data from    %
% West al. (2018) using an optimization based upon Approximate Bayesian   %
% Computation and described in West et al. (2019). Many of these analyses %
% are based upon scripts included in publicly available toolboxes and     %
% generously provided by their respective authors. Please see             %
% 'ABC_dependencies' for their licenses.                                  %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Timothy West, Nuffield Department of Clinical Neurosciences, University %
% of Oxford; Wellcome Centre for Human Neuroscience, University College   %
% London. 2018-2020   

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
clear; close all

% Initialise Analysis
%organise file paths, and add dependencies to path. Initializes
%configuration structure 'R' used throughout the scripts
R = basalGangliaRhythms_AddPaths(); 
% Add configurations/settings to R
R = setupBasalGangliaModel(R);

%% Figure (1): Model fit and comparison with data
simlength = 256; % Set simulation length
modID = 10; % This is the index of the winning model
BAA_plotModelFit(R,modID,simlength);
plotDataComparison % plots the time series in figure 1A and C

%% Figure (2): Simulating Inactivation of Connectivity
simlength = 32;
fresh = 1; % This flag will resimulate if 'fresh', else will load precomputed data and plot
BAA_sim_lesionExp(R,modID,simlength,fresh)
BAA_sim_lesionExp_signChange(R,modID,simlength,fresh)
BAA_sim_lesionExp_cortical_rhythms(R,modID,simlength,fresh)

%% Perform sweep over connections and simulate
BAA_sim_ConnectionSweep_v2(R,modID,100,2)
%%
BAA_sim_ConnectionSweep_dualpaths(R,modID,100,2)
plotSweepSpectraWrapper_dualpaths(R)


%% Figure (3): Sweep over connections and plot spectra
R = plotSweepSpectraWrapper(R); % This doesnt just plot but also works out the vector Krange for spread of viable connections
plotSweepSpectraWrapper_M2_SI(R); % Plot Cortical Spectra

%% Figure (4): Analysis of competing loops
[R,m,permMod,xsimMod] = getSimModelData_Draw(R,modID,32,0,500);
save('permModtmp_inflated','permMod')
load('permModtmp_inflated','permMod')
computeLoopAnalysis(R,permMod)

%% Perfom burst simulations and analysis
% Simulates the burst data using the range of induced STN beta 10% to 190%
simulateBurstData(R); % gets'_bKF' data, scaled to beta 10,100,190%
R = computeBurstWrapper_V2(R);

%% Figure (5): Compute timelocked analyses- burst coincidence
BB.struccmap = linspecer(4);
fresh = 1;
burstCoincidenceCheck(R,fresh)
burstPropertiesCheck(R,fresh)
burstLockCheck(R,fresh)
%% Figure (6): Compute timelocked analyses- relative burst timings
OnsetEvolveAnalysisMaster(R)

%% Figure (6 and 7): Closed loop stimulation- M2 and STN simulation 
simlength = 80; %80;
fresh = 1;
BAA_sim_fakeCloseLoop_StateDependency_v2(R,modID,simlength,fresh);

%% Figure (8): Plot State dependency of ARCs
ClosedLoop_StateDependence_PLot






%% SCRIPT GRAVE %%%
% % This is a analysis of Beta inputs STN/M2
% simlength = 32;
% fresh = 0;
% [R,MP] = BAA_sim_betaInputs(R,modID,simlength,fresh);
% 
% % This is pertubation analysis to get PRCs- at the moment its not very
% % clear what were looking for with it
% % [R] = BAA_sim_PRC(R,MP,500,0);
% [R] = BAA_sim_BetaPropagation(R,160,1); %remember reference of Boba and Hamacher 2015 for ZTE
% BB.range.RP = linspace(-pi,pi,7);
% BB = computeBetaBurstRPStats(R,BB);
%
% % TimeLockedBetaPLVAnalysis(R,BB,xsimMod,AS)
% % 
% % %% Bifurc
% % BAA_sim_InputSweep_v2(R,modID,35,4)

% % simlength = 50;
% % fresh = 1;
% % modID = 10;
% % [R] = BAA_sim_fakeCloseLoop(R,modID,simlength,fresh);
