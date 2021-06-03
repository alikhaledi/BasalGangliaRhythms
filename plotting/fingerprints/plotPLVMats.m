function plotPLVMats(R,PLVMat,CON,cmap,ha)
sblot = [2 1 3;0 0 0;5 7 9];
 
for state = 2:3
%     if state == 2
%         band = 1;
%     elseif state == 3
%         band = 1;
%     end
     pA = indepTTest(PLVMat{1,state,CON}(:,:,1),PLVMat{1,1,CON}(:,:,1),PLVMat{1,state,CON}(:,:,2),PLVMat{1,1,CON}(:,:,2),PLVMat{1,state,CON}(:,:,3),PLVMat{1,1,CON}(:,:,3));
     pA(isnan(pA)) = 0;
     pB = indepTTest(PLVMat{2,state,CON}(:,:,1),PLVMat{2,1,CON}(:,:,1),PLVMat{2,state,CON}(:,:,2),PLVMat{2,1,CON}(:,:,2),PLVMat{2,state,CON}(:,:,3),PLVMat{2,1,CON}(:,:,3));
     pB(isnan(pB)) = 0;
     p = pA + pB;
     
    A = PLVMat{1,state,CON}(:,:,1)./PLVMat{1,1,CON}(:,:,1);
    A(isnan(A)) = 0;
    A(logical(eye(size(A)))) = nan;
    B = PLVMat{2,state,CON}(:,:,1)./PLVMat{2,1,CON}(:,:,1);
    B(isnan(B)) = 0;
    B(logical(eye(size(B)))) = nan;    
    AB = A+B;
    
    AB = AB.*(p<0.05); 
    AB(AB==0) = nan;
    
    if state == 3
        C = PLVMat{1,1,CON}(:,:,1);
        C(isnan(C)) = 0;
        C(logical(eye(size(C)))) = nan;
        
        D = PLVMat{2,1,CON}(:,:,1);
        D(isnan(D)) = 0;
        D(logical(eye(size(D)))) = nan;
        CD = C + D;
        %             subplot(3,3,5); %(state-1)+(CON-1))
        axes(ha(5));
        cmapmid = brewermap(256,'Greys');
        plotMats(CD,cmapmid,R);
        title('Base Model')
                caxis([0.25 1])

    end
    
    %             subplot(3,3,sblot(CON,state)); %(state-1)+(CON-1))
    axes(ha(sblot(CON,state)));
    
    plotMats(AB,cmap,R)
    title(R.statename{CON,state})
    
end


function a = plotMats(AB,cmap,R)
imagesc(AB,'AlphaData',~isnan(AB))
colormap(gca,cmap)
a = gca;
set(a, 'ydir', 'reverse');
c = colorbar;
ylabel(c, 'Change in Burst PLV')

        caxis([0.25 1.5])
a.XTick = 1:6;
a.XTickLabel =     R.chsim_name;
a.XTickLabelRotation = -45;
% a.XLabel.String = 'Sensing';
a.YTick = 1:6;
a.YTickLabel =     R.chsim_name;
% a.YLabel.String = 'Overlapping';
hold on
plot([7 0],[7 0],'k','LineWidth',2)
axis square;
grid on

