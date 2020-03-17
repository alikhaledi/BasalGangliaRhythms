function burstAnalyAddPaths()
switch getenv('computername')
    case 'SFLAP-2'
        usname = 'Tim';
        gitpath = ['C:\Users\' usname '\Documents\Work\GIT'];
        madpath = 'MATLAB_ADDONS';
        spmpath = 'C:\Users\Tim\Documents\spm12';
    case 'FREE'
        usname = 'twest';
        gitpath = '\Documents\Work\GitHub';
        madpath = 'Work\MATLAB ADDONS';
        spmpath = 'C:\spm12';
    case 'DESKTOP-94CEG1L'
        usname = 'timot';
        gitpath =  '\Documents\GitHub';
        madpath = 'C:\Users\timot\Documents\Work\MATLAB ADDONS';
        spmpath = 'C:\Users\timot\Documents\GitHub\spm12';
    case 'TIM_PC'
        usname = 'Tim';
        gitpath =  'D:\GITHUB';
        madpath = 'D:\MATLAB ADDONS';
        spmpath = 'D:\GITHUB\spm12-master';
end

pathCell = regexp(path, pathsep, 'split'); onPath = any(strcmpi(spmpath, pathCell));
if ~onPath; addpath(spmpath); spm eeg; close all; end
addpath([madpath '\ParforProgMon'])
addpath([madpath '\aboxplot'])
addpath([madpath '\ATvDFA-package'])
addpath([madpath '\bplot\'])
addpath([madpath '\Circular_Statistics_Toolbox'])
addpath([madpath '\DrosteEffect-BrewerMap-221b913'])
addpath([madpath '\export_fig'])
addpath([madpath '\FMINSEARCHBND'])
addpath([madpath '\linspecer'])
addpath([madpath '\MEG_STN_Project'])
addpath([madpath '\Neurospec\neurospec21'])
addpath(genpath([madpath '\ParforProgMon']))
addpath([madpath '\sigstar-master'])
addpath([madpath '\sort_nat'])
addpath([madpath '\SplitVec'])
addpath([madpath '\TWtools'])
addpath([madpath '\violin'])
addpath(genpath([spmpath '\toolbox\xjview96\xjview']))
addpath(genpath([gitpath '\SimAnneal_NeuroModel\sim_machinery']))
addpath(genpath([gitpath '\SimAnneal_NeuroModel\Projects\Rat_NPD\routine\rat_InDirect_ModelComp']))
addpath(genpath([madpath '\boundedline-pkg']))
addpath(genpath([gitpath '\BrewerMap']))
addpath(genpath([gitpath '\BurstToolbox']))
addpath(genpath([gitpath '\SimAnneal_NeuroModel\Projects\Rat_NPD\data']));
addpath(genpath([gitpath '\SimAnneal_NeuroModel\Projects\Rat_NPD\model_fx']));
addpath(genpath([gitpath '\SimAnneal_NeuroModel\Projects\Rat_NPD\priors']));