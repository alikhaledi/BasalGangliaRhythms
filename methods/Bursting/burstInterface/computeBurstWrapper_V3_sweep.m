function [R,BB] = computeBurstWrapper_V3_sweep(R)
rootan = [R.rootn 'data\' R.out.oldtag '\ConnectionSweep'];
close all;
for CON = 3; %[1 3]
    BB = [];
    load([rootan '\BB_' R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_xsim_F1.mat'],'xsim')
    load([rootan '\BB_' R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_ck_1_F1.mat'],'ck_1')
    for state = 1:31; %size(ck_1,2)-1
        %% Find Burst Epochs
        % Get Data
        X = xsim{state}{1};
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
            
            
            % compute burst PLV
            [statBurstPLVl{band,state,CON},~,] = burstPLV(XPhi,burstSelInds,band);
           
            % compute burst overlaps
            [statBurstOverl{band,state,CON} sampBurstOverl{band,state,CON}] = burstOverlap(XF,fsamp,frqz,minper,band);
            
            % Compute Connectivity/Spectral Matrix
            [Pw pHz] = pwelch(X',fsamp,[],fsamp,fsamp);
            specMat{band,state,CON} = Pw(pHz>4 & pHz<=48,:);
                
            connectMat{band,state,CON} = computePLVConMat(XPhi,band);
            disp([band state CON])
        end
    end
end

save([rootan '\BB_' R.out.tag '_stateSweepConnectivityMatrix.mat'],'specMat','connectMat','statBurstOverl','statBurstPLVl')


a = 1;
% end