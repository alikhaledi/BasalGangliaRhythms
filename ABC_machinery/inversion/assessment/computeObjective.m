function [ACC R2] = computeObjective(R,r2)
if ~isfield(R.SimAn,'scoreweight')
    R.SimAn.scoreweight = [1 1];
    warning('Combined score weightings not available, setting to Unity')
end
R2 = R.SimAn.scoreweight(1)*(r2);
ACC = (R.SimAn.scoreweight(1)*(r2)) - (R.SimAn.scoreweight(2)*R.Mfit.DKL);