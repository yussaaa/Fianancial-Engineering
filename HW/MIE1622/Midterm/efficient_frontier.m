function[frontier] = efficient_frontier(Na, mu, Q, rf, w_init, w_Sharpe, w_ERC)
% Compute Markowitz efficient frontier

% Add path to CPLEX
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio128\cplex\matlab\x64_win64');

% Random data for 10 stocks
n = Na;

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
% fprintf ('Minimum variance portfolio:\n');
% fprintf ('Solution status = %s\n', cplex1.Solution.statusstring);
% fprintf ('Solution value = %f\n', cplex1.Solution.objval);
% fprintf ('Return = %f\n', sqrt(ret_minVar));
% fprintf ('Standard deviation = %f\n\n', sqrt(var_minVar));

% Compute maximum return portfolio
cplex2 = Cplex('max_Return');
cplex2.Model.sense = 'maximize';
cplex2.addCols(mu, [], lb, ub);
cplex2.addRows(b, A, b);
cplex2.Param.lpmethod.Cur = 6; % concurrent algorithm
cplex2.Param.barrier.crossover.Cur = 1; % enable crossover
cplex2.DisplayFunc = []; % disable output to screen
cplex2.solve();

% Display maximum return portfolio
w_maxRet = cplex2.Solution.x;
var_maxRet = w_maxRet' * Q * w_maxRet;
ret_maxRet = mu' * w_maxRet;
% fprintf ('Maximum return portfolio:\n');
% fprintf ('Solution status = %s\n', cplex2.Solution.statusstring);
% fprintf ('Solution value = %f\n', cplex2.Solution.objval);
% fprintf ('Return = %f\n', sqrt(ret_maxRet));
% fprintf ('Standard deviation = %f\n\n', sqrt(var_maxRet));

% Target returns
targetRet = linspace(ret_minVar,ret_maxRet,20);

% Compute efficient frontier
cplex3 = Cplex('Efficient_Frontier');
cplex3.addCols(zeros(n,1), [], lb, ub);
cplex3.addRows(targetRet(1), mu', inf);
cplex3.addRows(b, A, b);
cplex3.Model.Q = 2*Q;
cplex3.Param.qpmethod.Cur = 6; % concurrent algorithm
cplex3.Param.barrier.crossover.Cur = 1; % enable crossover
cplex3.DisplayFunc = []; % disable output to screen

w_front = [];
for i=1:length(targetRet)
    cplex3.Model.lhs(1) = targetRet(i);
    cplex3.solve();
    w_front = [w_front cplex3.Solution.x];
    var_front(i) = w_front(:,i)' * Q * w_front(:,i);
    ret_front(i) = mu' * w_front(:,i);
end

% Equally weighted Port
w_eq = repmat(1/n,1,n);
% Equally weighted portfolio return and variance
ret_equal_weight = dot(mu, w_eq);
var_equal_weight = w_EQW' * Q * w_eq;

%Initial Portfolio 
ret_init_portf = dot(mu, w_init);
var_init_portf = w_init' * Q * w_init;

%Max Sharpe Ratio
ret_Sharpe = dot(mu, w_Sharpe);
var_Sharpe = w_Sharpe' * Q * w_Sharpe;

%Risk Free Asset
ret_rf = rf/252 ;
var_rf = 0;

%Eqully risk contribution
ret_ERC = dot(mu, w_ERC);
var_ERC = w_ERC' * Q * w_ERC;

%Levergeded equaly weight 
w_leqr = 2*w_ERC;
ret_LERC = 2 * ret_ERC - ret_rf;
var_LERC = w_leqr' * Q * w_leqr;
end

rnd_portf = zeros(1000,n);

for i =1:1000
    rnd_portf(j,:) = rand(1,20);
end





% % Plot efficient frontier
% figure(1);
% set(gcf, 'color', 'white');
% plot(sqrt(var_front), ret_front, 'k-', 'LineWidth', 3)
% hold on;
% plot(sqrt(var_minVar), ret_minVar, 'rd', 'MarkerSize', 6)
% hold on;
% plot(sqrt(var_maxRet), ret_maxRet, 'ms', 'MarkerSize', 6)
% hold on;
% plot(sqrt(diag(Q)), mu, 'b.', 'MarkerSize', 18)
% xlabel('Standard deviation');
% ylabel('Expected return');
% title('Efficient Frontier')
% legend('efficient frontier', 'minimum variance portfolio', 'maximum return portfolio', 'individual stocks', 'Location', 'SouthWest')