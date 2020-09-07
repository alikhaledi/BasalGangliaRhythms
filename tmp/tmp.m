tData=6;  % in Secs
SR=500;
dtTime=1/SR;
nData=SR*tData;
tmaxis=(0:nData-1)*dtTime;

fr1=20;
am1=1;
am2=2;

dt1=am1*sin(2*pi*fr1*tmaxis+0);
dt2=am2*(rand(1,nData)-0.5);
dta=dt1+dt2;

%%
fhz=fr1;
BandWidth=2;

wo = fhz/(SR/2);  
bw = BandWidth/(SR/2);
qf =3;  % Q factor ( -qf);

[bnc,anc] = iirnotch(wo,bw,qf);
[bpk,apk] = iirpeak(wo,bw,qf);

%
%dt_nc=filtfilt(bnc,anc,dta);
%dt_pk=filtfilt(bpk,apk,dta);

dt_nc=filter(bnc,anc,dta);
dt_pk=filter(bpk,apk,dta);

% dt_pk=IIRPeak_Flt(dta,SR,fhz,BandWidth,qf);

lm1=min(dta); lm2=max(dta);
ux=figure; 
set(ux,'Position',[39 378 1477 618]);
subplot(3,1,1); plot(tmaxis,dta,'k');  ylim([lm1 lm2]); grid on; title('Native');
subplot(3,1,2); plot(tmaxis,dt_nc,'r'); ylim([lm1 lm2]); grid on; title('IIR Notch');
subplot(3,1,3); hold on;
plot(tmaxis,dt_pk,'g'); 
plot(tmaxis,(dta-dt_nc),'r'); 
ylim([lm1 lm2]); 
grid on; title('IIR Peak (green) and (Native - IIR Notch) red');
