function Rout = dataUpdateFix(Rorg,modID)
Rout = Rorg;
Rorg.out.dag = sprintf([Rorg.out.tag '_M%.0f'],modID);
[Rload,m,p] = loadABCData_160620(Rorg);
Rout.IntP.intFx = @ABC_fx_compile_120319; % update fx name
Rout.frqz = Rload.frqz;
Rout.data = Rload.data;