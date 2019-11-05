clc;
clear all;
format long

% Pricing a European option using Black-Scholes formula and Monte Carlo simulations
% Pricing a Barrier option using Monte Carlo simulations

S0 = 100;     % spot price of the underlying stock today
K = 105;      % strike at expiry
mu = 0.05;    % expected return
sigma = 0.2;  % volatility
r = 0.05;     % risk-free rate
T = 1.0;      % years to expiry
Sb = 110;     % barrier


% Implement your Black-Scholes pricing formula
[call_BS_European_Price, putBS_European_Price] = BS_european_price(S0, K, T, r, sigma);

% Implement your multi-step Monte Carlo pricing procedure for European
% option with different number of steps 
for  numPaths = 100:800:4000
for (numSteps = 1:25:1000)      % step number list
    i = (numSteps +25-1)/25;
    [callMC_European_Price_multi_step, putMC_European_Price_multi_step] = MC_european_price(S0, K, T, r, mu, sigma, numSteps, numPaths);
    step(i,1)=numSteps;
    diff(i,1) = callMC_European_Price_multi_step - call_BS_European_Price + putMC_European_Price_multi_step - putBS_European_Price;
    
end

% find index
diff = abs(diff);
index = find(diff == min(diff(:,1)));
disp(['Number of Steps =', num2str(numPaths)])
disp(['Min Difference between Monte Carlo and Black-Scholes pricing formula happened at step ',num2str(step(index,1))])
disp(['Difference between Monte Carlo and Black-Scholes pricing formula is ',num2str(diff(index,1))])
end