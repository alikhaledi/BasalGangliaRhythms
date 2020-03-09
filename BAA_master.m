%% Master - Compute Beta Burst Aligned Analysis for Indirect Comp %%%
clear; close all
R = basalGangliaRhythms_AddPaths();

%% Setup
%% Set Routine Pars
R = setupBasalGangliaModel(R);
modID = 10;
%% Plot Model Fit
simlength = 128;
BAA_plotModelFit(R,modID,simlength);

%% Inactivation Experiments
% BAA_sim_lesionExp_cortical_rhythms(R,modID,32,1)
simlength = 32;
fresh = 1;
BAA_sim_lesionExp(R,modID,simlength,fresh)
fresh = 1;
BAA_sim_lesionExp_cortical_rhythms(R,modID,simlength,fresh)

%% Perform Loop Analysis
[R,m,permMod,xsimMod] = getSimModelData_Draw(R,modID,32,0,500);
% This Draw is not working!
save('permModtmp_inflated','permMod')
load('permModtmp_inflated','permMod')

computeModelDrawParameters(R,permMod)


%% Plot Model Sweep Spectra
BAA_sim_ConnectionSweep_v2(R,modID,100,2)
R = plotSweepSpectraWrapper(R); % This doesnt just plot but also works out the vectors for spread
plotSweepSpectraWrapper_M2_SI(R); % Plot Cortical Spectra

%% Compute the Propagation Timing
% Make Clear which sets of data are being used.
simulateBurstData(R); % basically just gets'_bKF' data, scaled to beta 10,100,190%
R = computeBurstWrapper(R);
% Plot TimeLocked Statistics
BB.struccmap = linspecer(4);
% TimeLockAnalysisMaster(R); %,BB,[15 17 20]); % [15 17 20] for STN_GPE [1 6 8]
OnsetEvolveAnalysisMaster(R)

%% Compute Closed Loop Experiment
simlength = 50;
fresh = 0;
modID = 10;
[R] = BAA_sim_fakeCloseLoop(R,modID,simlength,fresh);

simlength = 70;
fresh = 1;
[R] = BAA_sim_fakeCloseLoop_StateDependency(R,modID,simlength,fresh);



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