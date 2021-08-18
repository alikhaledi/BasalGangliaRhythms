main = 1/18; % period (s)
driftFrq = 20;
drift = 1/driftFrq; % period (s)

error = drift-main; % (s)

errorCyc = 2*pi*driftFrq*error; % rad
errorDeg = rad2deg(errorCyc)

