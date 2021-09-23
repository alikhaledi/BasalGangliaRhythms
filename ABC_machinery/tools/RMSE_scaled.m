function [MSE RMSE] = RMSE_scaled(y,yhat)
E = (y - yhat);    % Errors
SE = E.^2; %(y - yhat).^2   % Squared Error

Epr = (y - mean(yhat));
SEpr = Epr.^2;
MSE = mean(SE)/mean(SEpr); %mean((y - yhat).^2)   % Mean Squared Error
RMSE = sqrt(MSE); %sqrt(mean((y - yhat).^2));  % Root Mean Squared Error
