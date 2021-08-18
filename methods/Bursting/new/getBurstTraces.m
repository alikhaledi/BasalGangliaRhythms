function [twin,m2Env,stnEnv,dPhi,sw_twin,sw_PLV,aftdPhi,aftdPhi2]= getBurstTraces(burstSelInds,winsize,fsamp,XEnv,XPhi)


for seg = 1:numel(burstSelInds)
    zeroP = burstSelInds{seg}(1); % This is burst Onset
    Win = zeroP-winsize(1):zeroP+winsize(2);
    aftWin = zeroP:zeroP+(winsize(2)/2);
    m2Env(:,seg) = XEnv(1,Win);
    stnEnv(:,seg) = XEnv(4,Win);
    

    
    dPhi(:,seg) = diff(XPhi([1 4],Win)); % M2/STN Phase
    sw_PLV(:,seg) = swPLV(dPhi(:,seg),0.15*fsamp); % Sliding Window PLV
    
    aftdPhi(:,seg) = diff(XPhi([1 4],aftWin));
    
    aftdPhi2(:,seg) = diff(XPhi([3 4],aftWin));
end
twin = 1000.*linspace(-winsize(1)/fsamp,winsize(2)/fsamp,size(m2Env,1));
sw_twin = 1000.*linspace(-winsize(1)/fsamp,winsize(2)/fsamp,size(sw_PLV,1));