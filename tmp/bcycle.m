clear; close all
load burstCycle

catInds = [burstSelInds{:}];
CTX = XF(1,catInds);
GPE = XF(3,catInds);
STN = XF(4,catInds);
 [N,Xedges,Yedges] = histcounts2(CTX,STN,50,'Normalization','pdf');
 midX = Xedges(1:end-1) + diff(Xedges);
  midY = Xedges(1:end-1) + diff(Yedges);

%  contour(midX,midY,N); hold on;
 
 
%  [~,sidx] = sort(cellfun(@numel,burstSelInds))
 [~,sidx] =sort(cellfun(@(x) mean(XEnv(4,x)),burstSelInds,'UniformOutput',1))

X1Traj = nan(1301,numel(burstSelInds)); % nan(max(cellfun(@numel,burstSelInds)),numel(burstSelInds));
X2Traj = nan(1301,numel(burstSelInds)); % nan(max(cellfun(@numel,burstSelInds)),numel(burstSelInds));
YTraj = nan(1301,numel(burstSelInds)); %nan(max(cellfun(@numel,burstSelInds)),numel(burstSelInds));
ATraj = nan(1301,numel(burstSelInds)); %nan(max(cellfun(@numel,burstSelInds)),numel(burstSelInds));

for i = 1:numel(burstSelInds)
    X = XF(1,burstSelInds{i});
            zci = find(diff(sign(X))>1,1,'first'); % location of last positive zero crossing
    zstep = burstSelInds{i}(1)+zci;
    zind = zstep-500:zstep+1200;
    
    X1 = XF(1,zind);
    X1Traj(1:numel(X1),i) = X1;
    X2 = XF(3,zind);
    X2Traj(1:numel(X2),i) = X2;
   
    Y = XF(4,zind);
    YTraj(1:numel(Y),i) = Y;

    A = XEnv(4,zind);
    ATraj(1:numel(A),i) = A;
end
tvec = -500:1:1200;
% X = reshape(XTraj(1:1288,1),[],14);
% Y = reshape(YTraj(1:1288,1),[],14);
% re

% cellfun(@(a) [a(1) a(end)],burstSelInds,'UniformOutput',0)
% plot(XTraj,YTraj)

subplot(2,1,1)
imagesc(tvec,1:100,(YTraj(:,sidx))')
ylabel('Burst # (ordered short to long')
xlabel('Time to onset (ms) (realigned)')
title('M2')

subplot(2,1,2)
imagesc(tvec,1:100,(X2Traj(:,sidx))')
ylabel('Burst # (ordered short to long')
xlabel('Time to onset (ms) (realigned)')
title('STN')

figure
subplot(2,1,1)
plot3(X2Traj(:,sidx(50)),YTraj(:,sidx(50)),ATraj(:,sidx(50)))
xlabel('STN'); ylabel('GPe'); zlabel('STN Beta Amplitude')
subplot(2,1,2)
plot3(X1Traj(:,sidx(50)),YTraj(:,sidx(50)),ATraj(:,sidx(50)))
xlabel('STN'); ylabel('M2'); zlabel('STN Beta Amplitude')

imagesc(sin(angle(hilbert(X1Traj)))-sin(angle(hilbert(YTraj))))


figure
subplot(2,1,1)
plot3(tvec,sin(angle(hilbert(X2Traj(:,sidx(1))))),sin(angle(hilbert(YTraj(:,sidx(1))))))
subplot(2,1,2)
plot3(tvec,sin(angle(hilbert(X1Traj(:,sidx(1))))),sin(angle(hilbert(YTraj(:,sidx(1))))))




figure
plot(sin(angle(hilbert(X2Traj(:,sidx(1))))),sin(angle(hilbert(YTraj(:,sidx(1))))))
