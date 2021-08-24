function [R,BB] = computeBurstWrapper_V3(R)
rootan = [R.rootn 'data\' R.out.tag '\ConnectionSweep'];
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
Xpair = [1 4];
for CON = 1:2
    BB = [];
    load([rootan '\BB_' R.out.tag '_DiscreteData.mat'],'dataSelect','dataProperties')
    
    for state = 1:4
        %% Find Burst Epochs
        % Get Data
        if state<4
            X = dataSelect{CON}{state}{1};
        else
            X  = dataSelect{CON}{1}{1};
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
    
    % figure(300+CON)
    % bplotBurstLength(segA,segL,conname{CON},statecmap{CON},CON)
    % set(gcf,'Position',[680   159   901   819])

    figure(431);
    subplot(2,1,splist(CON))
    cmap = statecmap{CON};
    for band = 1
        for state = [2 3]
            X = derv{band,state};
            Y = log10(intAmp{band,state}); %-log10(mean(intAmp{band,state}));
            scatter(X,Y,120,'MarkerFaceColor',cmap(state,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.3); hold on
            [r,p] = corrcoef(X,Y)
            if p(2)<0.05
                [xCalc yCalc b Rsq] = linregress(X',Y')
                plot(xCalc,yCalc,'LineWidth',2,'Color',cmap(state,:));
            end
        end
    end
    xlabel('Average relative phase instability'); ylabel('Burst power'); grid on; box off
    axis square
    set(gcf,'Position',[658    52   314   756])
    
        figure(CON)
        plotBurstTraces(twin,m2Env,stnEnv,sw_twin,sw_PLV,aftdPhi,middPhii,befdPhi,dPhi,statecmap{CON})
        set(gcf,'Position',[24 42 1813 756])
    
    figure(100)
    if CON == 1
        cmap = brewermap(128,'RdBu');
    elseif CON == 2
        cmap = brewermap(128,'*PRGn');
    end
    
    figure(101)
    plotOverlapMats(R,statBurstOverl,CON,cmap,ha)
    
    figure(102)
    %     plotConnectivityMats(R,connectMat,diffConnectCI,CON,cmap)
    plotConnectivityMats(R,statBurstPLVl,diffConnectCI,CON,cmap)
    % plotPLVMats(R,statBurstPLVl,CON,cmap,hb)
    
end
figure(101)
set(gcf,'Position',[401         321        1368         657])
figure(102)
set(gcf,'Position',[401         321        1368         657])

figure(300)
plotOverlapBars(sampBurstOverL)
save([rootan '\BB_' R.out.tag '_stateConnectivityMatrix.mat'],'specMat','connectMat','statBurstOverl','statBurstPLVl','segA','segL')

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
