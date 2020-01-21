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

% Lesion analysis
% BAA_sim_lesionExp_cortical_rhythms(R,modID,32,1)
simlength = 32;
fresh = 1;
BAA_sim_lesionExp(R,modID,simlength,fresh)

%% Plot Model Sweep Spectra
BAA_sim_ConnectionSweep_v2(R,modID,100,1)
R = plotSweepSpectraWrapper(R); % You need the R produced here!
plotSweepSpectraWrapper_M2_SI(R); % Plot Cortical Spectra

% BAA_sim_ConnectionSweep_v2(R,modID,500,2) % Compute smaller list

%% ISSUES IN THESE
% Make Clear which sets of data are being used.
simulateBurstData(R);
R = computeBurstWrapper(R);
%% Plot TimeLocked Statistics
BB.struccmap = linspecer(4);
TimeLockAnalysisMaster(R); %,BB,[15 17 20]); % [15 17 20] for STN_GPE [1 6 8]
OnsetEvolveAnalysisMaster(R)
%% 

% This is a analysis of Beta inputs STN/M2
[R,MP] = BAA_sim_betaInputs(R,10,32);

% This is pertubation analysis to get PRCs- at the moment its not very
% clear what were looking for with it
% [R] = BAA_sim_PRC(R,MP,500,0);
[R] = BAA_sim_fakeCloseLoop(R,500,1);

[R] = BAA_sim_BetaPropagation(R,160,1); %remember reference of Boba and Hamacher 2015 for ZTE
 
 
% BB.range.RP = linspace(-pi,pi,7);
% BB = computeBetaBurstRPStats(R,BB);
% 
% % TimeLockedBetaPLVAnalysis(R,BB,xsimMod,AS)




