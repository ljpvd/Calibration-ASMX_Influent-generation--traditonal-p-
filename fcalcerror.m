function [res,m,mae,mrse,obs,pred,n] = fcalcerror(t_exp,y_exp,t,y)

[a b] = size(y_exp);       % 2D matrix, a = number of samples, b = number of compounds. 
[c e] = size(y);           % 3D matrix, c = number of sim samples, d = number of compounds, e = number of simulations.

% First bucle, select the simulated data based on the observed data time frequency.
    j=0;
    k=1;
    for j=0:1/24:a
        index = (find (t==j));
        y_pred(k,:) = [y(index,:)]
        t_pred(k,:) = [t(index)];
        k=k+1;
    end


% Second bucle, calculate the different quality error indexes.
n = index;
time_obs = t_exp(:,1);
time_pred = t_pred;
obs = y_exp(:,1);
pred = [t_pred, y_pred];

for h=1:e
% Residuals
res(:,h) =  pred(:,h)-obs;
% Mean of residuals
m(h) = sum(res(:,h))./n;
% Mean absolute error
mae(h) = sum(abs(res(:,h)))./n;
% Root mean squared error
mrse(h) = sqrt(sum(res(:,h).^2)/n);
end

end