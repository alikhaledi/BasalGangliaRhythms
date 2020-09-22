function [R,BB] = computeBurstWrapper_V3(R)
rootan = [R.rootn 'data\' R.out.oldtag '\ConnectionSweep'];
close all;
conname = {'HD','','PS'};
statename = {'Fitted','A','B'};

for CON = [1 3]
    BB = [];
    load([rootan '\BB_' R.out.tag '_DiscreteData.mat'],'dataSelect','dataProperties')
    
    for state = 1:3
        %% Find Burst Epochs
        % Get Data
        X = dataSelect{CON}{state}{1};
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
            
            [burstSelInds,XF,XEnv,XPhi] = simpleBurstDefine(X,fsamp,frqz,minper);
            
            % Compute Connectivity Matrix
            cmblist = combnk(1:size(XPhi,1),2);
            conmat = zeros(6);
            for i=1:size(cmblist,1)
                XL = diff(XPhi([cmblist(i,1),cmblist(i,2)],:)); % relative phase
                
                if band == 1 % rotate for each band so you can fit in one matrix
                    conmat(cmblist(i,1),cmblist(i,2)) = PLV(XL);
                else
                    conmat(cmblist(i,2),cmblist(i,1)) = PLV(XL);
                end
            end
            connectMat{band,state} = conmat;
            
            % Setup window in which to look
            winsize(1) = 0.3*fsamp;
            winsize(2) = 1*fsamp;
            burstSelInds = cropBurstSelection(burstSelInds,winsize,size(X,2));% Remove bursts at ends
            
            [twin{band,state},m2Env{band,state},stnEnv{band,state},dPhi{band,state},sw_twin{band,state},sw_PLV{band,state},aftdPhi{band,state},aftdPhi2{band,state}] = getBurstTraces(burstSelInds,winsize,fsamp,XEnv,XPhi);
            
        end
    end
    figure(CON)
    plotBurstTraces(twin,m2Env,stnEnv,sw_twin,sw_PLV,aftdPhi)
    set(gcf,'Position',[680         367        1036         611])
    
    figure(100+CON)
    for state = 2:3
        A = connectMat{1,state}./connectMat{1,1};
        A(isnan(A)) = 0;
        B = connectMat{2,state}./connectMat{2,1};
        B(isnan(B)) = 0;
        AB = A+B;
        AB(logical(eye(size(AB)))) = nan;
        
        subplot(1,2,state-1)
        imagesc(AB,'AlphaData',~isnan(AB))
        a = gca;
        set(a, 'ydir', 'reverse');
        c = colorbar;
        ylabel(c, 'Change in PLV')
        
        caxis([0.5 1.5])
        colormap(brewermap(128,'RdBu'))
        a.XTick = 1:6;
        a.XTickLabel =     R.chsim_name;
        a.XTickLabelRotation = -45;
        
        a.YTick = 1:6;
        a.YTickLabel =     R.chsim_name;
        a.XTickLabelRotation = -45;
        title(statename{state})
    end
end
a = 1;
% end