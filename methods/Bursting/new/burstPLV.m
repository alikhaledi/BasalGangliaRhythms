function [statBurstPLVL,sampBurstPLVL,diffCI] = burstPLV(XPhi,burstSelInds,band)

for seg = 1:numel(burstSelInds)
    for chI = 1:6
        for chJ = 1:6
            XDPhi = diff(XPhi([chI,chJ],burstSelInds{seg}));
            burstPLVl(seg,chI,chJ) = PLV(XDPhi);
        end
    end
end

% Do some statistics
winsize = fix(min(cellfun(@numel,burstSelInds)));
bootPLV = [];
for cond = 1:2
    for i = 1:500
        chi = randi(size(XPhi,2)-winsize);
        chj = randi(size(XPhi,2)-winsize);
        
        XPhi_i = XPhi(1,chi:chi+winsize);
        XPhi_j = XPhi(1,chj+winsize:-1:chj);
        
        bootPLV(i,cond) = abs(PLV(XPhi_i-XPhi_j));
    end
end

bootdiff = diff(bootPLV,1,2);
diffCI = [prctile(bootdiff,5) prctile(bootdiff,95)];

statBurstPLVL = nan(6,6,4);
if band == 1
    statBurstPLVL = triu(squeeze(nanmean(burstPLVl,1)));
%     statBurstPLVL(:,:,1) = triu(squeeze(nanmean(abs(burstPLVl),1)));
%     statBurstPLVL(:,:,2) = triu(squeeze(nanstd(abs(burstPLVl),[],1)));
%     statBurstPLVL(:,:,3) = triu(squeeze(sum(~isnan(burstPLVl))));
%     statBurstPLVL(:,:,4) = triu(squeeze(nanmean(burstPLVl,1)));
    sampBurstPLVL = (burstPLVl);
elseif band == 2
    statBurstPLVL = tril(squeeze(nanmean(burstPLVl,1)));
%     statBurstPLVL(:,:,1) = tril(squeeze(nanmean(abs(burstPLVl),1)));
%     statBurstPLVL(:,:,2) = tril(squeeze(nanstd(abs(burstPLVl),[],1)));
%     statBurstPLVL(:,:,3) = tril(squeeze(sum(~isnan(burstPLVl))));
%     statBurstPLVL(:,:,4) = tril(squeeze(nanmean(burstPLVl,1)));
    sampBurstPLVL = (burstPLVl);
end

