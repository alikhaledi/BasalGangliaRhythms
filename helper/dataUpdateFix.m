function [Rout,m,p] = dataUpdateFix(Rorg,modID)
Rout = Rorg;
Rorg.out.dag = sprintf([Rorg.out.tag '_M%.0f'],modID);
[Rload,m,p] = loadABCData_160620(Rorg);
Rout = Rload;
Rout.path = Rorg.path;
Rout.out = Rorg.out;
Rout.IntP.intFx = @ABC_fx_compile_120319; % update fx name THIS IS NEW: ABC_fx_compile_120319
    Rout.obs.outstates(1) = 1; % change to middle layer
    m.outstates{1} = [1 0 0 0 0 0 0 0];

data_for = Rout.data.feat_emp{1};
data_rev = [];
%% flip i->j
for i = 1:size(data_for,2)
    for j = 1:size(data_for,3)
        if i~=j
            data_rev(1,i,j,:,:) = data_for(1,j,i,:,:);
            data_for(1,i,j,:,:) = data_for(1,i,j,:,:);
        end
    end
end
Rout.data.feat_emp{1} = data_rev;


