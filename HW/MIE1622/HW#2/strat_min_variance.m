function  [x_optimal cash_optimal] = strat_min_variance(x_init, cash_init, mu, Q, cur_prices, r_rf, init_value)
% this stragety is trying to find minimiumm varianvce profilo

% Add path to CPLEX
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio128\cplex\matlab\x64_win64');

n= length(x_init);
% Optimization problem data
lb = zeros(n,1);
ub = inf*ones(n,1);
A  = ones(1,n);
b  = 1;

% Compute minimum variance portfolio
cplex1 = Cplex('min_Variance');
cplex1.addCols(zeros(n,1), [], lb, ub);
cplex1.addRows(b, A, b);
cplex1.Model.Q = 2*Q;
cplex1.Param.qpmethod.Cur = 6; % concurrent algorithm
cplex1.Param.barrier.crossover.Cur = 1; % enable crossover
cplex1.DisplayFunc = []; % disable output to screen
cplex1.solve();

% Display minimum variance portfolio
w_minVar = cplex1.Solution.x;
var_minVar = w_minVar' * Q * w_minVar;
ret_minVar = mu' * w_minVar;

% total capital at the last period
% portfolio=portfolio price*share+cash
capital_minivar = cur_prices * x_init + cash_init;
share = floor(w_minVar.*capital_minivar./cur_prices'); % share on each stock and round this number to integer

x_optimal = share; % assign the weight to output
transaction = cur_prices * abs((x_optimal - x_init)) * 0.005;
cash_optimal = capital_minivar - cur_prices * share - transaction;

    % check whether cash account has negative balance
    if cash_optimal < 0
        fprintf('      ALERT: Min_variance strategy - CASH ACCOUNT BELOW ZERO \n');
    end

end



