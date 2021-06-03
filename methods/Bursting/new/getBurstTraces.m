function [twin,m2Env,stnEnv,dPhi,sw_twin,sw_PLV,aftdPhi,befdPhi,middPhii,allEnv,derv,intAmp]= getBurstTraces(burstSelInds,winsize,winloc,fsamp,XEnv,XPhi,Xpair)
twin = 1000.*linspace(-0.5,0.75,2501);

for seg = 1:numel(burstSelInds)
    zeroP = burstSelInds{seg}(1); % This is burst Onset
    Win = zeroP-(0.5*fsamp):zeroP+(0.75*fsamp);
    derivWin = zeroP:zeroP+(0.5*fsamp);
    
    %     befWin = zeroP-winsize(1):zeroP;
    %     midWin = zeroP:zeroP+(winsize(2)/2);
    %     aftWin = zeroP+(winsize(2)/2):zeroP+(winsize(2));
    
    befWin = zeroP+winloc(1)-winsize(1):...
        zeroP+winloc(1)+winsize(1);
    midWin = zeroP+winloc(2)-winsize(2):...
        zeroP+winloc(2)+winsize(2);
    aftWin = zeroP+winloc(3)-winsize(3):...
        zeroP+winloc(3)+winsize(3);
    m2Env(:,seg) = XEnv(Xpair(1),Win);
    stnEnv(:,seg) = XEnv(Xpair(2),Win);
    
    allEnv(:,:,seg) = XEnv(:,Win);
    % Make a local window here at maximum effect?
    [a maxInd] = max(stnEnv(:,seg));
    maxP = Win(maxInd);
    midWin = maxP-winsize(2):maxP+winsize(2);
    
    dPhi(:,seg) = wrapTo2Pi(diff(XPhi(Xpair,Win))); % M2/STN Phase
    %         WinT = linspace(-500,750,numel(Win));
    %         plot(WinT,dPhi(:,seg));
    sw_PLV(:,seg) = swPLV(dPhi(:,seg),0.15*fsamp); % Sliding Window PLV
    
    befdPhi(:,seg) =  wrapTo2Pi(diff(XPhi(Xpair,befWin)));
    middPhii(seg) = circ_mean(wrapTo2Pi(diff(XPhi(Xpair,burstSelInds{seg})))');
    aftdPhi(:,seg) =  wrapTo2Pi(diff(XPhi(Xpair,aftWin)));
    %     aftdPhi2(:,seg) = diff(XPhi([3 4],aftWin));
    Yx = diff(unwrap(diff(XPhi(Xpair,derivWin))));
    derv(seg) = mean(abs(Yx));
    intAmp(seg) = sum(XEnv(Xpair(2),derivWin))/(numel(derivWin)*fsamp); %power in au/s-1
    %     figure(523);clf; plot(twin,unwrap(dPhi(:,seg) )); yyaxis right; plot(twin,stnEnv(:,seg))
end
twin = 1000.*linspace(-0.5,0.75,size(m2Env,1));
sw_twin = 1000.*linspace(-0.5,0.75,size(sw_PLV,1));


%%% GRAVEYARD
% MULTIVARIATE PLV
% HTS: A Measure of Multivariate Phase Synchrony Using Hyperdimensional Geometry
%     Xpair2 = [1 4; 1 3; 6 3; 6 4; 6 5];
%     tmp_HTS = [];
%     for mset = 1:size(Xpair2,1)
%         dPhi = diff(XPhi(Xpair2(mset,:),Win)); % M2/STN Phase
%         tmp_HTS2(:,mset,seg) = swPLV(dPhi,0.15*fsamp); % Sliding Window PLV
%         tmp_HTS(:,mset) = swPLV(dPhi,0.15*fsamp); % Sliding Window PLV
%     end
%     sw_PLV(:,seg) = sqrt(mean(abs(tmp_HTS).^2,2));
