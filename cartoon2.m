     Fs = 1e3;
     t = (0:1000)'/Fs;
     x = sin(2*pi*[50 150 250].*t);
     x = sum(x,2) + 0.001*randn(size(t));
     bandpass(x,[100,200],Fs)