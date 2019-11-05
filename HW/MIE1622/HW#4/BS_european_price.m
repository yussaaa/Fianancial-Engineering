function [call_BS_European_Price, putBS_European_Price] = BS_european_price(S0, K, T, r, sigma)
    t = 0; % current tine stemp
    d1 = 1/(sigma * sqrt(T-t))*(log2(S0/K)+(r+0.5*sigma^2)*(T-t));
    d2 = d1 - sigma*sqrt(T-t);
% Black Scholes equation 
% Call Option Price Calculation 
    call_BS_European_Price = normcdf(d1)*S0 - normcdf(d2)*K*exp(-r*(T-t));
% Put Option price Calculation 
    putBS_European_Price = normcdf(-d2)*K*exp(-r*(T-t))-normcdf(-d1)*S0;

end