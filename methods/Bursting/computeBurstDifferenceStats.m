ctype = 1
load([R.rootn '\data\CloseLoop_stateDependency\CloseLoop_burstStats_save_' num2str(ctype) '.mat'],'durStore','ampStore','ppcStore','siStore')

[a pSup] = min(squeeze(burAmp(1,2,:,1)));
[a pAmp] = max(squeeze(burAmp(1,2,:,1)));

[h p ci stats] = ttest2(durStore{1,pSup},durStore{2,pSup});
[nanmean(durStore{1,pSup})-nanmean(durStore{2,pSup}) diff(ci)./2 stats.df stats.tstat p 0.05/5]

[h p ci stats] = ttest2(ampStore{1,pSup},ampStore{2,pSup});
[nanmean(ampStore{1,pSup})-nanmean(ampStore{2,pSup}) diff(ci)./2 stats.df stats.tstat p 0.05/5]

[h p ci stats] = ttest2(log(siStore{1,pSup}(5:end)),log(siStore{2,pSup}(5:end)));
[nanmean(siStore{1,pAmp})-nanmean(siStore{2,pAmp}) diff(ci)./2 stats.df stats.tstat p 0.05/5]

[h p ci stats] = ttest2(log(ppcStore{1,pSup}(3:end)),log(ppcStore{2,pSup}(3:end)));
[nanmean(ppcStore{1,pAmp})-nanmean(ppcStore{2,pAmp}) diff(ci)./2 stats.df stats.tstat p 0.05/5]
