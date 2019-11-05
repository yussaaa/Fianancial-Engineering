function  [x_optimal cash_optimal] = strat_max_Sharpe(x_init, cash_init, mu, Q, cur_prices)
% Stragety on finding the portfolio with the highest Sharpe ratio 
% Sharpe ratio which is also known as reward-to-variability ratio. It
% measures the ecxess return(risk premium) per unit of risk(deviation) in
% an investment assest or a portfolio 

    % Add the Cplex to the system path 
    addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio128\cplex\matlab\x64_win64\');
    
    n = length(x_init)+1; % Since variable kappa is intridced into the slover to s
    risk_fre = 0.025/252; % the rate of risk_free rate per year given in the main funtion
    mu = mu'-risk_fre; %( ui-rf)
    
    Q = [Q zeros(20,1) ; zeros(1,21)]; % Construct the right dimension for Q
    A = [mu, 0; ones(1,20),-1]; %two constrains
    b = [1;0];
    lb = zeros(n,1);
    ub = inf*ones(n,1);

    % Introducing the model 
    cplex3 = Cplex('max_Sharpe');
    cplex3.addCols(zeros(n,1), [], lb, ub);
    cplex3.addRows(b, A, b);
    
    % Specifying quadratic minimizing goal
    cplex3.Model.Q = 2*Q;
    cplex3.Param.qpmethod.Cur = 6; % concurrent algorithm 
    cplex3.Param.barrier.crossover.Cur = 1; % enable crossover 
    cplex3.DisplayFunc = []; % disable output to screen 
    cplex3.solve();

    solution = cplex3.Solution.x; %21*1
    
    % The output of the solver is actually y = w*k
    % to get the weight, w = y/k
    w_maxSharp = solution(1:20,:) / solution(21,:); %20*1
    
    capital_minivar = cur_prices * x_init + cash_init;
    share = floor(w_maxSharp.*capital_minivar./cur_prices'); % share on each stock and round this number to integer
    x_optimal = share; % assign the weight to output
    transaction = cur_prices * abs((x_optimal - x_init)) * 0.005;
    cash_optimal = capital_minivar - cur_prices * share - transaction;
end
