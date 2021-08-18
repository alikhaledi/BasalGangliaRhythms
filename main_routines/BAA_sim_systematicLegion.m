function [R] = BAA_sim_systematicLegion(Rorg,modID,simtime,fresh)
% Comopute simulations by sweeping across data
[R,m,permMod,xsimMod{1}] = getSimModelData_v3(Rorg,modID,simtime);
p = permMod{1}.par_rep{1};
R = setSimTime(R,simtime); % This sets the simulatiomn time by modify R structure
R.Bcond = -1; % No contrasting conditions

% Connection spec [ex/inh src tar]
AIJ = {[1 2 1],[1 4 1],[1 6 1],[1 3 4],[2 3 2],[2 5 2],[2 4 3],[1 5 4],[2 6 5],[1 1 6]};
Alist = 1:numel(AIJ);
Gint = 14;
Glist1 = numel(AIJ)+1:1:numel(AIJ)+Gint;
Glist2 = Glist1(end) + 1;
Glist6 = Glist2(end) + 1;
ilist = [Alist Glist1 Glist2 Glist6];

% Give all timeseries the same input
uc = innovate_timeseries(R,m);
uc{1} = uc{1}.*sqrt(R.IntP.dt);
R.obs.trans.norm = 0; % No normalization of spectra

% Simulate Base Model
Pbase = p;
[r2,~,feat_sim,dum,xsim_gl] = computeSimData120319(R,m,uc,Pbase,0); % Simulates the new model
featBase = feat_sim;
% Find resulting power statistics of the simulated data:
% STN beta
[powI_B_base(1),peakI_B_base(1),freqI_B_base(1)] = findSpectralStats(R.frqz,squeeze(feat_sim(1,4,4,1,:)),[14 30]);
[powI_B1_base(1),peakI_B1_base(1),freqI_B1_base(1)] = findSpectralStats(R.frqz,squeeze(feat_sim(1,4,4,1,:)),[14 21]);
[powI_B2_base(1),peakI_B2_base(1),freqI_B2_base(1)] = findSpectralStats(R.frqz,squeeze(feat_sim(1,4,4,1,:)),[21 30]);
% M2 Beta
[powI_m2_B_base(1),peakI_m2_B_base(1),freqI_m2_B_base(1)] = findSpectralStats(R.frqz,squeeze(feat_sim(1,1,1,:)),[14 30]);
[powI_m2_B1_base(1),peakI_m2_B1_base(1),freqI_m2_B1_base(1)] = findSpectralStats(R.frqz,squeeze(feat_sim(1,1,1,1,:)),[14 21]);
[powI_m2_B2_base(1),peakI_m2_B2_base(1),freqI_m2_B2_base(1)] = findSpectralStats(R.frqz,squeeze(feat_sim(1,1,1,1,:)),[21 30]);


if fresh == 1
    zi = 0;
    for i = ilist
        zi = zi +1;
        PI{i} = Pbase; % Load parameter base
        if any(intersect(i,Alist))
            PI{i}.A{AIJ{i}(1)}(AIJ{i}(2),AIJ{i}(3)) = -32; % Remove the ith connection
        elseif any(intersect(i,Glist1))
            [dum ind] = intersect(Glist1,i);
            PI{i}.int{1}.G(ind)= -32;
        elseif any(intersect(i,Glist2))
            [dum ind] = intersect(Glist2,i);
            PI{i}.int{2}.G(ind)= -32;            
        elseif any(intersect(i,Glist6))
            [dum ind] = intersect(Glist6,i);
            PI{i}.int{6}.G(ind)= -32;            
        end
    end
    
    
    parfor i = ilist
        [r2,~,feat_sim,dum,xsim_gl] = computeSimData120319(R,m,uc,PI{i},0); % Simulates the new model
        featSave{i} = feat_sim;
        % Find resulting power statistics of the simulated data (STN
        % beta
        [powI_B(i),peakI_B(i),freqI_B(i)] = findSpectralStats(R.frqz,squeeze(feat_sim(1,4,4,1,:)),[14 30]);
        [powI_B1(i),peakI_B1(i),freqI_B1(i)] = findSpectralStats(R.frqz,squeeze(feat_sim(1,4,4,1,:)),[14 21]);
        [powI_B2(i),peakI_B2(i),freqI_B2(i)] = findSpectralStats(R.frqz,squeeze(feat_sim(1,4,4,1,:)),[21 30]);
        
        [powI_m2_B(i),peakI_m2_B(i),freqI_m2_B(i)] = findSpectralStats(R.frqz,squeeze(feat_sim(1,1,1,1,:)),[14 30]);
        [powI_m2_B1(i),peakI_m2_B1(i),freqI_m2_B1(i)] = findSpectralStats(R.frqz,squeeze(feat_sim(1,1,1,1,:)),[14 21]);
        [powI_m2_B2(i),peakI_m2_B2(i),freqI_m2_B2(i)] = findSpectralStats(R.frqz,squeeze(feat_sim(1,1,1,1,:)),[21 30]);
        
        fitIJ(i) = r2;
        disp([i])
    end
    
    % STN Low Beta
    X = (powI_B1-powI_B1_base)./(powI_B1_base).*100;
    conlist{1} = find((X)<-25);
    % STN High Beta
    X = (powI_B2-powI_B2_base)./(powI_B2_base).*100;
    conlist{2} = find(X<-25);
    
    % M2 Low Beta
    X = (powI_m2_B1-powI_m2_B1_base)./(powI_B1_base).*100;
    conlist{3} = find(X<-25);
    % M2 High Beta
    X = (powI_m2_B2-powI_m2_B2_base)./(powI_m2_B2_base).*100;
    conlist{4} = find(X<-25);
    
end