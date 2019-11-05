function [callMC_Barrier_Knockin_Price_1_step, putMC_Barrier_Knockin_Price_1_step] = MC_barrier_knockin_price(S0, Sb, K, T, r, mu, sigma, numSteps, numPaths)
    
    % sb barrier of the option 2
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
    
    % This is a up option since sb>s0
    knockinpth = paths; 
%     knockin(1,:) = K;
    for iPath = 1:numPaths
        for iStep = 2:numSteps+1
            if any(paths(iStep,iPath) < Sb)
               knockinpth(iStep,iPath) = K;
            end
        end
    end
        
    % Calculate the return of any Put option in Random Walk
    PutPayoffT = max(K-knockinpth(numSteps+1,:),0);

    % Calculate the return of any Call option in Random Walk
    CallPayoffT = max(knockinpth(numSteps+1,:)-K,0);

    % Discount back to the time of option 
    putMC_Barrier_Knockin_Price_1_step = mean(PutPayoffT)*exp(-r*T);
    callMC_Barrier_Knockin_Price_1_step = mean(CallPayoffT)*exp(-r*T);
    
%     %Plot paths
%     figure;
%     
%     set(gcf, 'color', 'white');
%     plot(0:numSteps, knockinpth', 'Linewidth', 1);
%     title('Geometric Random Walk Paths - MC Barrier option price', 'FontWeight', 'bold');
%     xlabel('number of step')
%     ylabel('Asset price')

end