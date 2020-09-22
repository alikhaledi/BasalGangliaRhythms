function coeff = PLV(phi)
coeff = abs(sum(exp(1i*phi)))./numel(phi);
end