function R = prepareRatData_NoGauss_Animal_NPD(R,animal,state)

load([R.filepathn '\NPD_paper_RatNPD_150618.mat']);
NPDmat = fyA;
load([R.filepathn '\nsPow_paper_RatNPD_150618.mat']);
nsPowmat = fyA;
load([R.filepathn '\Frq_paper_RatNPD_150618.mat']);
meannpd_data = nan(1,6,6,4,311);
% condsel = [1 2];
R.condnames = {'OFF'};
R.Bcond = 1;
condsel = state;
chsel = [1 3 2 4]'; % M1 STR GPe STN

for C =1:numel(R.condnames)
    %     X = squeeze(fyA(:,:,:,:,2,:)); {i,j,dirc,cond,sub}
    for i = 1:size(chsel,1)
        for j = 1:size(chsel,1)
            if i==j
                F_data = fxA(:,1);
                F_model = R.frqz;
                % Log transform of group average
                Pxy = abs(log10(nsPowmat{chsel(i),state,animal}));
                % Remove the powerline frequency
                Pxy(F_data>48 & F_data<52) = [];
                F_data(F_data>48 & F_data<52) = [];
                % Put in the same frequency space as the models
                Pxy = interp1(F_data,Pxy,F_model);
                % Regress out the 1/f background
                [xCalc yCalc b Rsq] = linregress(log10(F_model)',Pxy');
                    Pxy = Pxy-fliplr(yCalc');
                    Pxy = 10.^Pxy;
                % Normalize and zero-base
                Pxy = (Pxy-mean(Pxy))./std(Pxy);
                Pxy = Pxy-min(Pxy);
                meannpd_data(C,i,j,1,:) = Pxy;
            else
                
                for k = 1:size(NPDmat,3)
                    F_data = fxA(:,1);
                    F_model = R.frqz;
                    Cxy = NPDmat{chsel(i),chsel(j),k,state,animal};
                    
                    % Remove the powerline frequency
                    Cxy(F_data>48 & F_data<52) = [];
                    F_data(F_data>48 & F_data<52) = [];
                    % Put in the same frequency space as the models
                    Cxy = interp1(F_data,Cxy,F_model);
                    meannpd_data(C,i,j,k,:) = Cxy;
                end
                
            end
        end
    end
end

% Set data as working version
R.data.feat_emp = meannpd_data;
% squeeze(meannpd_data(1,1,1,1,:))
R.data.feat_xscale = R.frqz;

