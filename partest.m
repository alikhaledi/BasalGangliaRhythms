
a = nan(1e3,100,20);
parfor p = 1:20
    x = rand(1e3,100);
    a(:,:,p) = x./4;
    disp(p)
end