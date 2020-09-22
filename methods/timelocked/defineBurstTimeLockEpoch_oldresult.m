function TL = defineBurstTimeLockEpoch2(BB,TL,cond)
% EpochTime
Bo = 0;
preBo = [Bo(1)+ floor((TL.periodT(1)/1e3)*BB.fsamp):Bo(1)-1]; %pre burst onset (indices)
postBo = [Bo(1): Bo(1) + floor((TL.periodT(2)/1e3)*BB.fsamp)]; % post burst onset (indices)
epochdef = [preBo(1):postBo(end)];
TL.epochT = linspace(TL.periodT(1),TL.periodT(2),size(epochdef,2));

postBo2 = [Bo(1): Bo(1) + floor((TL.periodT(3)/1e3)*BB.fsamp)]; % post burst onset (indices)
epochdef_off = [preBo(1):postBo2(end)];
TL.epochT_off = linspace(TL.periodT(1),TL.periodT(3),size(epochdef_off,2));

% LocalEps
localeps = BB.epsAmpfull;
segInds = BB.segInds{cond};
% Work first with lengths
clear BEpoch REpoch PLVpoch dRPdEpoch meanPLV maxAmp minAmp epsCross usedinds
for i = 1:numel(segInds)-1
    %% Set up vector of indices for each burst pre and burst onset at STN
    Bo = segInds{i}; % This is the indices for this burst segment
    preBo = [Bo(1)+ floor((TL.periodT(1)/1e3)*BB.fsamp):Bo(1)-1]; %pre burst onset
    postBo = [Bo(1): Bo(1) + floor((TL.periodT(2)/1e3)*BB.fsamp)]; % post burst onset
    epochdef = [preBo(1):postBo(end)];
    postBo2 = [Bo(1): Bo(1) + floor((TL.periodT(3)/1e3)*BB.fsamp)]; % post burst onset (indices)
    epochdef_off = [preBo(1):postBo2(end)];
    
    %     plotExampleBurst % This will plot time series in figure 6A
    
    %% Now look at all bursts across the network
    if preBo(1)>0 && postBo2(end)<size(BB.AEnv{cond},2)
        % Find onset Time aligned to beta onset
        A = BB.AEnv{cond}(:,epochdef); % Amplitude Envelope for epoch [-150 to +150]
        AH = BB.AEnv{cond}(:,epochdef).*hanning(numel(epochdef))'; % with Hanning
        AOFF = BB.AEnv{cond}(:,epochdef_off).*hanning(numel(epochdef_off))'; % Envelope for [-150 to +150]
        % Find Crossing times with respect to STN onset
        epsCross = []; maxCross = []; epsLast = [];
        for L = 1:size(BB.AEnv{cond},1)
            if any(A(L,:)>localeps(L)) % For finding maximums locally
                % Find local bursts
                betaBurstInds = SplitVec(find(AOFF(L,:)>localeps(L)),'consecutive'); % Split up data based upon the target threshold
                segL = cellfun('length',betaBurstInds); % compute burst lengths
                
                sinit = [];
                for ci = 1:numel(segL)
                    sinit(ci) = betaBurstInds{ci}(1); % get burst initiations
                end
%                 sinit(segL<BB.period) = inf; % remove small bursts
%                 sinit(abs(sinit(:)-(size(preBo,2)+1))> floor((TL.periodT(2)/1e3)*BB.fsamp)) = inf; % remove bursts over initial window away
                if any(~isinf(sinit))
                    [dum closeBurst] = min(abs(sinit(:)-(size(preBo,2)+1))); % find burst closest to STN onset
                    ci = betaBurstInds{closeBurst}([1 end]);
%                     epsCross(L) =  TL.epochT_off(ci(1));
%                     epsLast(L) =  TL.epochT_off(ci(2));
                    ci(1) = betaBurstInds{1}(1);
                     ci(2) = betaBurstInds{end}(end);
                    epsCross(L) =  TL.epochT_off(ci(1));
                    epsLast(L) =  TL.epochT_off(ci(2));
                    
                else
                    epsCross(L) = NaN;
                    epsLast(L) = NaN;
                end
            else
                epsCross(L) = NaN;
                epsLast(L) = NaN;
                
            end
        end
        %         TL.maxT{cond}(:,i) = maxCross;
        TL.onsetT{cond}(:,i) = epsCross;
        TL.onsetOffT{cond}(:,i) = epsLast;
        TL.segDur{cond}(:,i) = epsLast-epsCross;
        A = (A-min(A,2))./std(A,[],2);
        %         A = A.^2; %.^2;
        TL.amp{cond}(:,:,i) = A;
        raw = BB.data{cond}(:,epochdef); %.*hanning(numel(epochdef))'; % amplitude data
        TL.raw{cond}(:,:,i) =raw;
        BP = BB.BP{cond}(:,epochdef); %.*hanning(numel(epochdef))'; % amplitude data
        TL.BP{cond}(:,:,i) = BP;
        
        % STR_M2 RP
        phi = []; dPhi = [];
        for L = 1:size(BB.Phi{cond},1)
            phi(L,:) = BB.Phi{cond}(L,epochdef);
            uphi(L,:) = unwrap(phi(L,:));
            dPhi(L,:) = [NaN abs((diff(uphi(L,:))))];
        end
        TL.phi{cond}(:,:,i) = phi;
        TL.dPhi{cond}(:,:,i) = dPhi;
        
        
        
    end
end


%% SCRIPT GRAVE
% Get Indices of Long Bursts
%         TL.durFlag{cond}(i,1) = BB.segDur{cond}(i) > prctile(BB.segDur{1},75);
%         TL.durFlag{cond}(i,2) = BB.segDur{cond}(i) < prctile(BB.segDur{1},25);
%         uwRP = unwrap(BB.Phi{cond}(1,:))-unwrap(BB.Phi{cond}(2,:));
%         RP = wrapToPi(uwRP(epochdef));
%         TL.STR_M2_RP{cond}(i) = circ_mean(RP');
%         TL.STR_M2_dRP{cond}(:,i) = diff(uwRP(epochdef));

%         PLVbase = nanmedian(BB.PLV{cond});
%         [dum T(1)] = min(abs(BB.SWTvec{cond}-BB.TSw(epochdef(1))));
%         T(2) = T(1) + floor(sum(abs(periodT/1000))/diff(BB.TSw(1:2)));
%         if epochdef(end)<size(BB.AEnv{cond},2) && epochdef(1) > 0 && T(2)<=size(BB.PLV{cond},2)
%             BEpoch(:,:,i) = 1*zscore(BB.AEnv{cond}(:,epochdef),0,2).*hanning(numel(epochdef))'; % ch x time x burstN
%             REpoch(:,:,i) = 1*zscore(BB.Tvec{cond}(:,epochdef),0,2).*hanning(numel(epochdef))';
% %             BEpoch(:,:,i) = 0.2*BB.A{cond}(:,epochdef); %.*hanning(numel(epochdef))'; % ch x time x burstN
% %             REpoch(:,:,i) = 0.5*BB.rawTime{cond}(:,epochdef).*hanning(numel(epochdef))';
%             PLVpoch(:,i) = 100*(BB.PLV{cond}(1,T(1):T(2))-PLVbase)/PLVbase ;
%             dRPdEpoch(:,i) = dRPdt(epochdef)';
%             meanPLV(i) = mean(PLVpoch(:,i)); %computePPC(squeeze(BB.Phi([1 4],Bo)));
%             maxAmp(i) = max(BB.AEnv{cond}(4,Bo));
%             minAmp(i) = min(BB.AEnv{cond}(4,preBo));
%             Segpoch(i) = segL(segInds(i));
%         end