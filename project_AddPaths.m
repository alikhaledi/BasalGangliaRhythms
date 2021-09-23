function R = project_AddPaths(R)
%
switch getenv('computername')
    case 'DESKTOP-1QJTIMO'
        gitpath =  'C:\Users\Tim West\Documents\GitHub';
        madpath = 'C:\Users\Tim West\Documents\MATLAB ADDONS';
         spmpath = 'C:\Users\Tim West\Documents\GitHub\spm12';
    case 'DESKTOP-0HO6J14'
        gitpath =  'D:\GITHUB';
        madpath = 'C:\Users\timot\OneDrive\Documents\Work\MATLAB ADDONS';
%         gitpath2 = 'D:\GITHUB_STORE';
    % case 'PC-NAME'
    % gitpath = 'path to general github repo
    % madpath = 'path to dependencies'
    otherwise
        error('You need to add your PC specific paths - ...please look at project_AddPaths.m for template')
end
% Sometimes the main project is some other folder to Github default
if ~exist('gitpath2','var')
    gitpath2 = gitpath;
end

R.rootn = [gitpath2 '\BasalGangliaRhythms\'];

% Add the root
addpath(genpath(R.rootn))

% Add SPM/Fieldtrip
pathCell = regexp(path, pathsep, 'split'); onPath = any(strcmpi(spmpath, pathCell));
if ~onPath; addpath(spmpath); spm eeg; close all; end

% Add External Toolboxes
addpath(genpath([gitpath2 '\ABC_Inference_Neural_Paper\sim_machinery'])); % load manually for now
% addpath(genpath([gitpath '\BurstToolbox']))
addpath([madpath '\Circular_Statistics_Toolbox'])
% addpath([madpath '\Neurospec\neurospec21'])
addpath([madpath '\DrosteEffect-BrewerMap-221b913'])
addpath([madpath '\SplitVec'])
addpath([madpath '\permutest'])
addpath([madpath '\tight_subplot'])
addpath([madpath '\TWtools'])


