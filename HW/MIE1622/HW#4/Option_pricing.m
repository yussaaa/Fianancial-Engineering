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

numSteps =  12;
numPaths = 1000;

% Define variable numSteps to be the number of steps for multi-step MC
% numPaths - number of sample paths used in simulations

%% Implement your Black-Scholes pricing formula
[call_BS_European_Price, putBS_European_Price] = BS_european_price(S0, K, T, r, sigma);

%% Implement your one-step Monte Carlo pricing procedure for European option
[callMC_European_Price_1_step, putMC_European_Price_1_step] = MC_european_price(S0, K, T, r, mu, sigma, 1, numPaths);

% Implement your multi-step Monte Carlo pricing procedure for European option
[callMC_European_Price_multi_step, putMC_European_Price_multi_step] = MC_european_price(S0, K, T, r, mu, sigma, numSteps, numPaths);

%% Implement your one-step Monte Carlo pricing procedure for Barrier option
% numSteps = 1;
[callMC_Barrier_Knockin_Price_1_step, putMC_Barrier_Knockin_Price_1_step] = ...
    MC_barrier_knockin_price(S0, Sb, K, T, r, mu, sigma, numSteps, numPaths);


% Implement your multi-step Monte Carlo pricing procedure for Barrier option
[callMC_Barrier_Knockin_Price_multi_step, putMC_Barrier_Knockin_Price_multi_step] = ...
    MC_barrier_knockin_price(S0, Sb, K, T, r, mu, sigma, numSteps, numPaths);
%%
disp(['Black-Scholes price of an European call option is ',num2str(call_BS_European_Price)])
disp(['Black-Scholes price of an European put option is ',num2str(putBS_European_Price)])
disp(['One-step MC price of an European call option is ',num2str(callMC_European_Price_1_step)])
disp(['One-step MC price of an European put option is ',num2str(putMC_European_Price_1_step)])
disp(['Multi-step MC price of an European call option is ',num2str(callMC_European_Price_multi_step)])
disp(['Multi-step MC price of an European put option is ',num2str(putMC_European_Price_multi_step)])
disp(['One-step MC price of an Barrier call option is ',num2str(callMC_Barrier_Knockin_Price_1_step)])
disp(['One-step MC price of an Barrier put option is ',num2str(putMC_Barrier_Knockin_Price_1_step)])
disp(['Multi-step MC price of an Barrier call option is ',num2str(callMC_Barrier_Knockin_Price_multi_step)])
disp(['Multi-step MC price of an Barrier put option is ',num2str(putMC_Barrier_Knockin_Price_multi_step)])
%% Difference between MC and MC barrier
disp(' ')
disp('--- one step ---')
disp(['Difference between European call option and Barrier call option is ',num2str(callMC_European_Price_1_step - callMC_Barrier_Knockin_Price_1_step)])
disp(['Difference between European put option and Barrier put option is  ',num2str(putMC_European_Price_1_step - putMC_Barrier_Knockin_Price_1_step)])
disp(' ')
disp('--- multi-step ---')
disp(['Difference between European call option and Barrier call option is ',num2str(callMC_European_Price_multi_step - callMC_Barrier_Knockin_Price_multi_step)])
disp(['Difference between European put option and Barrier put option is  ',num2str(putMC_European_Price_multi_step - putMC_Barrier_Knockin_Price_multi_step)])

%% barrier option with different sigma
% Implement your one-step Monte Carlo pricing procedure for Barrier option
disp(' ')
disp('--- sigma = 0.1 ---')
sigma = 0.1;  % volatility decreased 10%
% one-step calculation
[callMC_Barrier_Knockin_Price_1_step, putMC_Barrier_Knockin_Price_1_step] = ...
    MC_barrier_knockin_price(S0, Sb, K, T, r, mu, sigma, 1, numPaths);

% Implement your multi-step Monte Carlo pricing procedure for Barrier option
[callMC_Barrier_Knockin_Price_multi_step, putMC_Barrier_Knockin_Price_multi_step] = ...
    MC_barrier_knockin_price(S0, Sb, K, T, r, mu, sigma, numSteps, numPaths);

disp('when volatility decreased by 10%: sigma = 0.1')
disp(['One-step MC price of an Barrier call option is ',num2str(callMC_Barrier_Knockin_Price_1_step)])
disp(['One-step MC price of an Barrier put option is ',num2str(putMC_Barrier_Knockin_Price_1_step)])
disp(['Multi-step MC price of an Barrier call option is ',num2str(callMC_Barrier_Knockin_Price_multi_step)])
disp(['Multi-step MC price of an Barrier put option is ',num2str(putMC_Barrier_Knockin_Price_multi_step)])

disp(' ')
disp('--- sigma = 0.3 ---')
sigma = 0.3;  % volatility increased 10%
[callMC_Barrier_Knockin_Price_1_step, putMC_Barrier_Knockin_Price_1_step] = ...
    MC_barrier_knockin_price(S0, Sb, K, T, r, mu, sigma, 1, numPaths);

% Implement your multi-step Monte Carlo pricing procedure for Barrier option
[callMC_Barrier_Knockin_Price_multi_step, putMC_Barrier_Knockin_Price_multi_step] = ...
    MC_barrier_knockin_price(S0, Sb, K, T, r, mu, sigma, numSteps, numPaths);

disp('when volatility increased by 10%: sigma = 0.3')
disp(['One-step MC price of an Barrier call option is ',num2str(callMC_Barrier_Knockin_Price_1_step)])
disp(['One-step MC price of an Barrier put option is ',num2str(putMC_Barrier_Knockin_Price_1_step)])
disp(['Multi-step MC price of an Barrier call option is ',num2str(callMC_Barrier_Knockin_Price_multi_step)])
disp(['Multi-step MC price of an Barrier put option is ',num2str(putMC_Barrier_Knockin_Price_multi_step)])

%% Plot results
% figure(1);
%%%%%%%%%%% Insert your code here %%%%%%%%%%%%
for (S0 = 0:10:250)
    i = (S0+10)/10;
    [call_BS_European_Price, putBS_European_Price] = BS_european_price(S0, K, T, r, sigma);
    [callMC_European_Price_multi_step, putMC_European_Price_multi_step] = MC_european_price(S0, K, T, r, mu, sigma, numSteps, numPaths);
    
    asset_price(i,1) = S0;
    call_BS(i,1)=call_BS_European_Price;
    put_BS(i,1)=putBS_European_Price;
    callMC_Eu(i,1)=callMC_European_Price_multi_step;
    putMC_Eu(i,1)=putMC_European_Price_multi_step;

end

% %% Plot price 
% figure;
% set(gcf, 'color', 'white');
% plot(asset_price, call_BS, 'Color', 'r','Linewidth', 2);
% hold on;
% plot(asset_price, callMC_Eu, 'Color', 'b','Linewidth', 2);
% title('Call option value', 'FontWeight', 'bold');
% legend('BS european price', 'MC european price (multi-step)');
% xlabel('Asset Price')
% ylabel('Option Value')
% hold off;
% 
% figure;
% set(gcf, 'color', 'white');
% plot(asset_price, put_BS,'Color', 'r', 'Linewidth', 2);
% hold on;
% plot(asset_price, putMC_Eu,'Color', 'b', 'Linewidth', 2);
% title('Put option value', 'FontWeight', 'bold');
% legend('BS european price', 'MC european price (multi-step)');
% xlabel('Asset Price')
% ylabel('Option Value')
% hold off;