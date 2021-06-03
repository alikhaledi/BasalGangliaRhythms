function bonfSigStat(X,base,fy,permnum,sigh,ofset,cmap)
ssList = 1:100:size(base,1);
p = [];
for sl = 1:numel(ssList)
    [h p(sl)] = ttest2(base(ssList(sl),:),fy(ssList(sl),:));
end

p(p>(0.01/numel(ssList))) = NaN;
a = sigh + (ofset);
sigY = a.*~isnan(p);
sigY(sigY==0) = NaN;
plot(X(ssList),sigY,'Color',cmap,'LineWidth',2)
