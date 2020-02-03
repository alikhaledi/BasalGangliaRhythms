function [bpowr_br,fpow_br,bpowr,fpow,bcohr,fcoh] = computeBetaSpectralStats(Hz,feat)
bpow = []; fpow = [];
for ck = 1:numel(feat)
    [dum b] = max(feat{ck}(1,4,4,3,Hz>14 & Hz<21)); % Low Beta Power
    fpow_br(ck) = Hz(b) + Hz(1);
    bpowr_br(ck) =  sum(feat{ck}(1,4,4,3,Hz>14 & Hz<21))./sum(Hz>14 & Hz<30);
    
    
    [dum b] = max(feat{ck}(1,4,4,3,Hz>14 & Hz<30)); % Full Beta Power
    fpow(ck) = Hz(b) + Hz(1);
    bpowr(ck) = sum(feat{ck}(1,4,4,3,Hz>14 & Hz<30))./sum(Hz>14 & Hz<30); % Full Beta Power
    
    
    [dum b] = max(feat{ck}(1,4,1,4,Hz>14 & Hz<30)); % Full Beta M2/STN Coh
    fcoh(ck) = Hz(b) + Hz(1);
    bcohr(ck) = mean(feat{ck}(1,4,1,4,Hz>14 & Hz<30));
end