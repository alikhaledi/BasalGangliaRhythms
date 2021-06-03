close all; clear
subplot(2,1,1)
fsamp = 1000;
dt = 1/fsamp;
frq = 14;
t = 0:dt:0.6;
aEnvMod = exp(-t/0.1)- exp(-t./0.01);
tPre = 0:dt:0.1;
tNew = [tPre tPre(end)+dt:dt:t(end)+tPre(end)+dt];
aEnvMod = [repmat(0.01,1,numel(tPre)) aEnvMod+0.01];

plot(tNew,aEnvMod)
hold on
% aEnvModStoc = aEnvMod + 0.1*randn(size(aEnvMod));


B = aEnvMod.*sin(2*pi*frq*tNew);
BStoc = B + 0.0*randn(size(aEnvMod));
plot(tNew,BStoc)

d = designfilt('bandpassfir','FilterOrder',120, ...
         'CutoffFrequency1',5,'CutoffFrequency2',50, ...
         'SampleRate',fsamp);
     

BFilt = filtfilt(d,BStoc);
plot(tNew,BFilt)

subplot(2,1,2)
sclz = -1;
Env = 1./(abs(hilbert(BFilt)).^2);
dscz = 1/(1/sum(Env));
for i = 1:numel(BFilt)-1
    sclz(i+1) = sclz(i) + (Env(i)/dscz);
end
sclz(end)

% sclz = linspace(-2,2,numel(BFilt)); %.*abs(hilbert(BFilt));
plot(sclz+imag(hilbert(BFilt)),real(hilbert(BFilt)))