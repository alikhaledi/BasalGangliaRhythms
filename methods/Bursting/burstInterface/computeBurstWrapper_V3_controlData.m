function [R,BB] = computeBurstWrapper_V3_controlData(R)
% rootan = [R.rootn 'data\ConnectionSweep'];
rootan = 'D:\Data\ratdata_050816';
close all;
% Overlap subplots
figure(101)
ha =   tight_subplot(3,3,0.05);
delete(ha([2 4 6 8]))
splist = [1 2]; %subplot list

scmap = brewermap(4,'Set1');
statecmap{1} = [0 0 0; scmap(1:2,:)];
statecmap{2} = [0 0 0; scmap(3:4,:)];

R.statename = {'Fitted','HD-Down','HD-Up'; 'Fitted','PS-Down','PS-Up'};
anim = {'L4','L6','L13','L18','L19','L20','L22','L23','L24'};
Xpair = [1 4];
CON =3; state = 1;
BB = [];
for aniN = 1:3
    load([rootan '\' anim{aniN} '.mat'],'dataSelect','dataProperties')
    
    state = 1;
    %% Find Burst Epochs
    % Get Data
    if state<4
        X = dataSelect{CON}{state}{1};
    else
        X  = dataSelect{CON}{1}{1}; % state 4 is used for permutation comparisons
    end
    for band = 1:2
        if band == 1
            fhz = 18; butf = [14 21];
        elseif band == 2
            fhz = 24; butf = [21 30];
        end
        
        % Set parameters
        fsamp = 1/R.IntP.dt;
        frqz = butf;
        minper = 3; % Min 3 cycles
        %         peakHz = dataProperties(4,1,CON); % This loads in the peak STN frequency
        
        [burstSelInds,XF,XEnv,XPhi,epsAmp,segL{band,state,CON},segA{band,state,CON}] = simpleBurstDefine(X,fsamp,frqz,minper);
        
        if state == 4
            burstSelInds = permuteInds(burstSelInds,X);
        end
        
        
        % compute burst PLV
        [statBurstPLVl{band,state,CON},~,] = burstPLV(XPhi,burstSelInds,band);
        
        % compute burst overlaps
        [statBurstOverl{band,state,CON},sampBurstOverL{band,state,CON}] = burstOverlap(XF,fsamp,frqz,minper,band);
        
        
        [connectMat{band,state,CON},diffConnectCI{band,state,CON}] = computePLVConMat(XPhi,band);
        
        % Compute Connectivity/Spectral Matrix
        [Pw pHz] = pwelch(X',fsamp,[],fsamp,fsamp);
        Pw = (Pw-mean(Pw(:)))./std(Pw(:));
        specMat{band,state,CON} = Pw(pHz>4 & pHz<=48,:);
        
        % Setup window in which to look
        winsize(1) = 0.075*fsamp;
        winsize(2) = 0.075*fsamp; %0.15*fsamp;
        winsize(3) = 0.075*fsamp;
        winsize(4) = 0.075*fsamp;
        winloc(1) =  -0.35*fsamp;
        winloc(2) =  0*fsamp;
        winloc(3) =  +0.5*fsamp;
        
        burstSelInds = cropBurstSelection(burstSelInds,winsize,size(X,2));% Remove bursts at ends
        
        [twin{band,state},m2Env{band,state},stnEnv{band,state},...
            dPhi{band,state},sw_twin{band,state},sw_PLV{band,state},...
            aftdPhi{band,state},befdPhi{band,state},middPhii{band,state},~,derv{band,state},intAmp{band,state}] = getBurstTraces(burstSelInds,winsize,winloc,fsamp,XEnv,XPhi,Xpair);
    end
end

function burstSelIndsPerm = permuteInds(burstSelInds,X)
L = 1:size(X,2); flag = 1;
for i = 1:numel(burstSelInds)
    flag = 1;
    while flag
        p = randi(numel(L));
        pinds = p:p+numel(burstSelInds{i});
        if (pinds(end)+2000)<size(X,2) &&  (pinds(1)-2000)>1
            flag = 0;
        end
    end
    burstSelIndsPerm{i} = pinds;
end
