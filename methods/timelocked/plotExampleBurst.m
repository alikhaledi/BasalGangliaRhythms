TL.periodT = [-1000 1000];
gap = 10;
cmap = linspecer(4);
localeps = prctile(abs(hilbert(BB.BP{1}))',75);

for  i = 4:numel(BB.segInds{1})
    Bo = segInds{i};
    preBo = [Bo(1)+ floor((TL.periodT(1)/1e3)*BB.fsamp):Bo(1)-1]; %pre burst onset
    postBo = [Bo(1): Bo(1) + floor((TL.periodT(2)/1e3)*BB.fsamp) + 1]; % post burst onset
    epochdef = [preBo(1):postBo(end)];
    TL.epochT = linspace(TL.periodT(1),TL.periodT(2),size(epochdef,2));
    
    p = plot(TL.epochT,BB.BP{1}(:,epochdef)'+ [0:-gap:-(gap*3)]);
    for l = 1:4
        p(l).Color = cmap(l,:);
    end
    hold on
    p = plot(TL.epochT,abs(hilbert(BB.BP{1}(:,epochdef)'))+ [0:-gap:-(gap*3)]);
    for l = 1:4
        p(l).Color = cmap(l,:);
    end
    
    p = plot(TL.epochT,repmat(localeps',1,numel(TL.epochT))' + [0:-gap:-(gap*3)]);
    for l = 1:4
        p(l).Color = cmap(l,:);
    end
    
    pause
    clf
end

plot([0 200],[-35 -35])
ylim([-40 5])

a = 1;