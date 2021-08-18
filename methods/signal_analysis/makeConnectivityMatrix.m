function conmat = makeConnectivityMatrix(XPhi,band)
cmblist = combnk(1:size(XPhi,1),2);
conmat = zeros(size(XPhi,1));
for i=1:size(cmblist,1)
    XL = diff(XPhi([cmblist(i,1),cmblist(i,2)],:)); % relative phase
    
    if band == 1 % rotate for each band so you can fit in one matrix
        conmat(cmblist(i,1),cmblist(i,2)) = PLV(XL);
    else
        conmat(cmblist(i,2),cmblist(i,1)) = PLV(XL);
    end
end