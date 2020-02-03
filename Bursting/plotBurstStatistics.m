close all

for CON = [1 3]
    R.condname{1} = '10%';
    R.condname{10} ='100%';
    R.condname{19} = '190%';
    rootan = [R.rootn 'data\' R.out.oldtag '\ConnectionSweep'];
    
    load([rootan '\BBA_' R.out.tag '_Sims_CON_' num2str(CON) '_bKF.mat'],'BB')
    
    
    R.condcmap = linspecer(30);
    %     R.condname([1 6 8 15 17 20]) =  {'1% M2->STN','Fitted','150% M2->STN','1% STN ->GPe','Fitted','150% STN->GPe'};
    R.condnames =  R.condname;
    
    % Setup Bin Ranges
    BB.range.Amp = linspace(0,50,20);% 1:0.25:5; %0:5:120; % Group: 0:3:120; single: 0:5:80
    BB.range.Amplr = linspace(0,50,10);% 1:0.25:5; %0:5:120; % Group: 0:3:120; single: 0:5:80
    % BB.range.Dur = linspace(50,1400,20);
    BB.range.Dur = linspace((20),(1800),20);
    BB.range.Durlr = linspace(log10(20),log10(1800),10);
    
    BB.range.segDur = linspace(1.5,3.2,24); %0:100:1800; % Group: 0:100:1800; single: 25:150:1800
    BB.range.AmpPrc = 0:5:100;
    
    % Setup Plot Limits
    BB.plot.lims.burfreq = [0 35; 0 35];
    BB.plot.lims.PLV = [-100 100]; %[0 0.35];
    BB.plot.lims.PLVun = [0 0.6];
    BB.plot.lims.wPLV = [-10 10];
    BB.plot.lims.Amp =  [0 30];
    BB.plot.lims.wAmpPrc =  [-2 6];
    BB.plot.lims.Dur = [0 1600]; %log10([15 1800]);
    BB.plot.durlogflag = 0;
    % Compute  Burst Amplitude/Duration Statistics
    % for i = 1:5; F(i) = figure; end
    %     if CON == 1
    condsel = [1 10 19]; % conditions to be selected
    R.condcmap(condsel,:) =  linspecer(3);
    
    
    figure
    BB = AmpDurStatistics(R,BB,condsel,CON);
end