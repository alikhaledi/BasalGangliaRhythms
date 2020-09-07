function R = basalGangliaRhythms_AddPaths(R)

switch getenv('computername')
    case 'SFLAP-2'
        usname = 'Tim'; gitpath = '\Documents\Work\GIT'; madpath = 'MATLAB_ADDONS';
        spmpath = 'C:\Users\Tim\Documents\spm12';
    case 'FREE'
        usname = 'twest'; gitpath = '\Documents\Work\GitHub'; madpath = 'Work\MATLAB ADDONS';
        spmpath = 'C:\spm12';
    case 'DESKTOP-94CEG1L'
        usname = 'timot';
        gitpath =  'C:\Users\timot\Documents\GitHub';
        gitpath2 =  'D:\GITHUB_STORE\';
        madpath = 'C:\Users\timot\Documents\Work\MATLAB ADDONS';
        spmpath = 'C:\Users\timot\Documents\GitHub\spm12';
    case 'TIM_PC'
        usname = 'Tim';
        gitpath =  'D:\GITHUB';
        madpath = 'D:\MATLAB ADDONS';
        spmpath = 'D:\GITHUB\spm12-master';
    case 'DESKTOP-1QJTIMO'
        usname = 'Tim West';
        gitpath =  'C:\Users\Tim West\Documents\GitHub';
        madpath = 'C:\Users\Tim West\Documents\MATLAB ADDONS';
        spmpath = 'C:\Users\Tim West\Documents\GitHub\spm12';
end
% Sometimes the main project is some other folder to Github default
if ~exist('gitpath2','var')
    gitpath2 = gitpath;
end

R.rootn = [gitpath2 '\BasalGangliaRhythms\'];

% Add the root
addpath(genpath(R.rootn))

% Add External Toolboxes
addpath(genpath([gitpath '\ABC_Inference_Neural_Paper\sim_machinery']))
pathCell = regexp(path, pathsep, 'split'); onPath = any(strcmpi(spmpath, pathCell));
if ~onPath; addpath(spmpath); spm eeg; close all; end
addpath(genpath([gitpath '\BurstToolbox']))
addpath([madpath '\Circular_Statistics_Toolbox'])
addpath([madpath '\Neurospec\neurospec21'])
addpath([madpath '\DrosteEffect-BrewerMap-221b913'])
addpath([madpath '\SplitVec'])
addpath([madpath '\TWtools'])
% pathCell = regexp(path, pathsep, 'split'); onPath = any(strcmpi(spmpath, pathCell));
% if ~onPath; addpath(spmpath); spm eeg; close all; end
% addpath(['C:\Users\' usname '\Documents\' madpath '\ParforProgMon'])
% addpath(['C:\Users\' usname '\Documents\' madpath '\aboxplot'])
% addpath(['C:\Users\' usname '\Documents\' madpath '\ATvDFA-package'])
% addpath(['C:\Users\' usname '\Documents\' madpath '\bplot\'])
% addpath(['C:\Users\' usname '\Documents\' madpath '\Circular_Statistics_Toolbox'])
% addpath(['C:\Users\' usname '\Documents\' madpath '\cirHeatmap'])
% addpath(['C:\Users\' usname '\Documents\' madpath '\DrosteEffect-BrewerMap-221b913'])
% addpath(['C:\Users\' usname '\Documents\' madpath '\export_fig'])
% addpath(['C:\Users\' usname '\Documents\' madpath '\FMINSEARCHBND'])
% addpath(['C:\Users\' usname '\Documents\' madpath '\HotellingT2'])
% addpath(['C:\Users\' usname '\Documents\' madpath '\linspecer'])
% addpath(['C:\Users\' usname '\Documents\' madpath '\MEG_STN_Project'])
% addpath(['C:\Users\' usname '\Documents\' madpath '\Neurospec\neurospec21'])
% addpath(genpath(['C:\Users\' usname '\Documents\' madpath '\ParforProgMon']))
% addpath(['C:\Users\' usname '\Documents\' madpath '\sigstar-master'])
% addpath(['C:\Users\' usname '\Documents\' madpath '\sort_nat'])
% addpath(['C:\Users\' usname '\Documents\' madpath '\SplitVec'])
% addpath(['C:\Users\' usname '\Documents\' madpath '\TWtools'])
% addpath(['C:\Users\' usname '\Documents\' madpath '\violin'])
% addpath(genpath(['C:\spm12\toolbox\xjview96\xjview']))
% addpath(genpath(['C:\Users\' usname '\' gitpath '\SimAnneal_NeuroModel\sim_machinery']))
% addpath(genpath(['C:\Users\' usname '\' gitpath '\SimAnneal_NeuroModel\Projects\Rat_NPD\routine\rat_InDirect_ModelComp']))
% addpath(genpath(['C:\Users\' usname '\Documents\' madpath '\boundedline-pkg']))
% addpath(genpath(['C:\Users\' usname '\' gitpath '\BrewerMap']))
% addpath(genpath(['C:\Users\' usname '\' gitpath '\BurstToolbox']))
% addpath(genpath(['C:\Users\' usname '\' gitpath '\highdim']))
% addpath(genpath(['C:\Users\' usname '\' gitpath '\SimAnneal_NeuroModel\Projects\Rat_NPD\data']));
% addpath(genpath(['C:\Users\' usname '\' gitpath '\SimAnneal_NeuroModel\Projects\Rat_NPD\model_fx']));
% addpath(genpath(['C:\Users\' usname '\' gitpath '\SimAnneal_NeuroModel\Projects\Rat_NPD\priors']));