function plotOverlapMats(R,overlapMat,CON,cmap,ha)
sblot = [2 1 3;0 0 0;5 4 6];
 
for state = 2:3
%     if state == 2
%         band = 1;
%     elseif state == 3
%         band = 1;
%     end
     pA = indepTTest(overlapMat{1,state,CON}(:,:,1),overlapMat{band,1,CON}(:,:,1),overlapMat{band,state,CON}(:,:,2),overlapMat{band,1,CON}(:,:,2),overlapMat{band,state,CON}(:,:,3),overlapMat{band,1,CON}(:,:,3));
     pA(isnan(pA)) = 0;
     
     p = pA;
     
    A = overlapMat{band,state,CON}(:,:,1)./overlapMat{band,1,CON}(:,:,1);
    A(isnan(A)) = 0;
    A(logical(eye(size(A)))) = nan;
    A = A.*(p<0.05); 
    A(A==0) = nan;
    
    if state == 3
        C = overlapMat{band,1,CON}(:,:,1);
        C(isnan(C)) = 0;
        C(logical(eye(size(C)))) = nan;
        %             subplot(3,3,5); %(state-1)+(CON-1))
        axes(ha(sblot(CON,1)));
        cmapmid = brewermap(256,'Greys');
        plotMats(C,cmapmid,R);
        title('Base Model')
                caxis([0.25 1])

    end
    
    %             subplot(3,3,sblot(CON,state)); %(state-1)+(CON-1))
    axes(ha(sblot(CON,state)));
    
    plotMats(A,cmap,R)
    title(R.statename{CON,state})
    
end


function plotMats(AB,cmap,R)
imagesc(AB,'AlphaData',~isnan(AB))
colormap(gca,cmap)
a = gca;
set(a, 'ydir', 'reverse');
c = colorbar;
ylabel(c, 'Change in Burst Overlap')

        caxis([0.25 1.5])
a.XTick = 1:6;
a.XTickLabel =     R.chsim_name;
a.XTickLabelRotation = -45;
a.XLabel.String = 'Sensing';
a.YTick = 1:6;
a.YTickLabel =     R.chsim_name;
a.YLabel.String = 'Overlapping';
% hold on
% plot([7 0],[7 0],'k','LineWidth',2)
axis square;
grid on

