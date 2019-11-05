function [callMC_European_Price_1_step, putMC_European_Price_1_step] = MC_european_price(S0, K, T, r, mu, sigma, numSteps, numPaths)
% Monte Carlo simulation - using Geometric Random Walk with constant drift and volatility
    
    % mu = annual drift, sigma = annual volatility
    % T is the total length of time for the path (in years)
    % dT is the time increment (in years)
       
    paths = zeros(numSteps+1, numPaths);
    dT = T/numSteps;
    
    % Vector of paths will store realizations of the asset price
    % First asset price is the initial price
    paths(1,:) = S0;
 
    % Generate paths
    for iPath = 1:numPaths
        for iStep = 1:numSteps
            paths(iStep+1, iPath) = paths(iStep, iPath) * exp( (mu - 0.5*sigma^2)*dT + sigma*sqrt(dT)*normrnd(0,1) );
        end
    end 

    % Calculate the return of any Put option in Random Walk
    PutPayoffT = max(K-paths(numSteps+1,:),0);

    % Calculate the return of any Call option in Random Walk
    CallPayoffT = max(paths(numSteps+1,:)-K,0);

    % Discount back to the time of option 
    putMC_European_Price_1_step = mean(PutPayoffT)*exp(-r*T);
    callMC_European_Price_1_step = mean(CallPayoffT)*exp(-r*T);
    
%     %Plot paths
%     figure;
%     set(gcf, 'color', 'white');
%     plot(0:numSteps, paths', 'Linewidth', 1);
%     title('Geometric Random Walk Paths - MC simulated price', 'FontWeight', 'bold');
%     xlabel('number of step')
%     ylabel('Asset price')

end