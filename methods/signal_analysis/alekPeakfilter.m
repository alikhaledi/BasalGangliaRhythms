function XF = alekPeakfilter(X,fhz,SR)
BandWidth=3;
wo = fhz/(SR/2);
bw = BandWidth/(SR/2);
qf =3;  % Q factor ( -qf);
% [bnt,ant] = iirnotch(wo,bw,qf);
            [bpk,apk] = iirpeak(wo,bw,qf);
%                         [dum,bpk,apk] = ft_preproc_bandpassfilter(X, fsamp, butf,4,'but','twopass');
%                         XF =filtfilt(bpk,apk,X')';

XF_nt =filtfilt(bnt,ant,X')';
XF = XF_nt;

end