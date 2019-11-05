clc;
clear all;
format long

% CSV file with price data
input_file_prices  = 'Daily_closing_prices.csv';

% Read daily prices
if(exist(input_file_prices,'file'))
  fprintf('\nReading daily prices datafile - %s\n', input_file_prices)
  fid = fopen(input_file_prices);
     % Read instrument tickers
     hheader  = textscan(fid, '%s', 1, 'delimiter', '\n');
     headers = textscan(char(hheader{:}), '%q', 'delimiter', ',');
     tickers = headers{1}(2:end);
     vheader = textscan(fid, '%[^,]%*[^\n]');
     dates = vheader{1}(1:end);
  fclose(fid);
  data_prices = dlmread(input_file_prices, ',', 1, 1);
else
  error('Daily prices datafile does not exist')
end

% Convert dates into array [year month day]
format_date = 'mm/dd/yyyy';
dates_array = datevec(dates, format_date);
     % Read time periods
dates_array = dates_array(:,1:3);

% Remove datapoints for year 2014
day_ind_start0 = 1;
day_ind_end0 = length(find(dates_array(:,1)==2014));
data_prices = data_prices(day_ind_end0+1:end,:);
dates_array = dates_array(day_ind_end0+1:end,:);
dates = dates(day_ind_end0+1:end,:);

% Compute means and covariances for Question 2
day_ind_start = 1;
day_ind_end = 39;
cur_returns = data_prices(day_ind_start+1:day_ind_end,:) ./ data_prices(day_ind_start:day_ind_end-1,:) - 1;
mu = mean(cur_returns)';  % Expected returns for Question 2
Q = cov(cur_returns);     % Covariances for Question 2


% Question 1

% Specify quantile level for VaR/CVaR
alf = 0.95;

% Positions in the portfolio
positions = [100 0 0 0 0 0 0 0 200 500 0 0 0 0 0 0 0 0 0 0]';

% Number of assets in universe
Na = size(data_prices,2);

% Number of historical scenarios
Ns = size(data_prices,1);

%%%%% Insert your code here 
loss_1d = - diff(data_prices);

for i = 1:Ns +1 -10
    loss_10d(i,:) = data_prices(i,:) - (data_prices(i+10 -1,:));   
end

% Number of losses in the given data in 1 day situation 
Ns_1 = Ns - 1;
% % Number of losses in the given data in 10 day situation 
Ns_10 = Ns - 10 +1;

%Portfolio loss per 1d in monetory value 
portf_loss_1d = sort(loss_1d*positions);
portf_loss_10d = sort(loss_10d*positions);

%Var Cvar calculation using historical data
%Var at quantile level=alf
VaR1  = portf_loss_1d(ceil(Ns_1 * alf));
VaR10 = portf_loss_10d(ceil(Ns_10 * alf));
%CVar at quantile level=alf
CVaR1 = (1/(Ns_1*(1-alf))) * ( (ceil(Ns_1*alf)-Ns_1*alf) * VaR1 + sum(portf_loss_1d(ceil(Ns_1*alf)+1:Ns_1)) );
CVaR10 = (1/(Ns_10*(1-alf))) * ( (ceil(Ns_10*alf)-Ns_10*alf)* VaR10 + sum(portf_loss_10d(ceil(Ns_10*alf)+1:Ns_10))); 


%Var Cvar calculation by using simulated normal distribution 
% VaR1_n = mean(portf_loss_1d) + norminv(alf,0,1)*std(portf_loss_1d);
% VaR10_n = mean(portf_loss_10d) + norminv(alf,0,1)*std(portf_loss_10d);
% 
% CVaR1_n = mean(portf_loss_1d) + (normpdf(norminv(alf,0,1))/(1-alf))*std(portf_loss_1d);
% CVaR10_n = mean(portf_loss_10d) + (normpdf(norminv(alf,0,1))/(1-alf))*std(portf_loss_10d);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mean and variance calculated by using mean and variance of -return rate
% which is the loss 
%%%%%% For daily loss 
day_ind_start = 1;
day_ind_end = 503;
cur_returns = data_prices(day_ind_start+1:day_ind_end,:) ./ data_prices(day_ind_start:day_ind_end-1,:) - 1;
mu_1 = -mean(cur_returns)';  
Q_1 = cov(cur_returns);   

value = data_prices(1,:) * positions;
w_p = (data_prices(1,:) .* positions')' / value;
mu_p = mu_1' * w_p * value;
std_p = sqrt(w_p' * Q_1 * w_p)  * value;

VaR1_n = mu_p + norminv(alf,0,1)*std_p;
CVaR1_n = mu_p + (normpdf(norminv(alf,0,1))/(1-alf))*std_p;

%%%%%% For loss per 10 days 
day_ind_start_10 = 1;
day_ind_end_10 = 495;
cur_returns_10 = data_prices(day_ind_start+9:day_ind_end,:) ./ data_prices(day_ind_start:day_ind_end-9,:) - 1;
mu_10 = -mean(cur_returns_10)';  
Q_10 = cov(cur_returns_10);   

value = data_prices(1,:) * positions;
w_p = (data_prices(1,:) .* positions')' / value;
mu_p_10 = mu_10' * w_p * value;
std_p_10 = sqrt(w_p' * Q_10 * w_p) * value;

VaR10_n = mu_p_10 + norminv(alf,0,1)*std_p_10;
CVaR10_n = mu_p_10 + (normpdf(norminv(alf,0,1))/(1-alf))*std_p_10;

fprintf('Historical 1-day VaR %4.1f%% = $%6.2f,   Historical 1-day CVaR %4.1f%% = $%6.2f\n', 100*alf, VaR1, 100*alf, CVaR1)
fprintf('    Normal 1-day VaR %4.1f%% = $%6.2f,       Normal 1-day CVaR %4.1f%% = $%6.2f\n', 100*alf, VaR1_n, 100*alf, CVaR1_n)
fprintf('Historical 10-day VaR %4.1f%% = $%6.2f,   Historical 10-day CVaR %4.1f%% = $%6.2f\n', 100*alf, VaR10, 100*alf, CVaR10)
fprintf('    Normal 10-day VaR %4.1f%% = $%6.2f,       Normal 10-day CVaR %4.1f%% = $%6.2f\n', 100*alf, VaR10_n, 100*alf, CVaR10_n)


% Plot a histogram of the distribution of losses in portfolio value for 1 day 
figure(1)
[frequencyCounts, binLocations] = hist(portf_loss_1d, 100);
bar(binLocations, frequencyCounts); hold on;
line([VaR1 VaR1], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
%normf = ( 1/(std(portf_loss_1d)*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(portf_loss_1d))/std(portf_loss_1d)).^2 ); normf = normf * sum(frequencyCounts)/sum(normf);
%plot(binLocations, normf, 'r', 'LineWidth', 3); hold on;
line([CVaR1 CVaR1], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-.');
hold off;
text(1.25*VaR1, max(frequencyCounts)/1.9, 'CVaR1'); text(0.7*VaR1_n, max(frequencyCounts)/1.9, 'VaR1');
title('Daily Loss horizon')
xlabel('Loss')
ylabel('Number of Scenario')
% Plot a histogram of the distribution of losses in portfolio value for 10 days
figure(2)
[frequencyCounts, binLocations] = hist(portf_loss_10d, 100);
bar(binLocations, frequencyCounts); hold on;
line([VaR10 VaR10], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
%normf = ( 1/(std(portf_loss_10d)*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(portf_loss_10d))/std(portf_loss_10d)).^2 ); normf = normf * sum(frequencyCounts)/sum(normf);
%plot(binLocations, normf, 'r', 'LineWidth', 3); hold on;
line([CVaR10 CVaR10], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-.');
hold off;
text(0.98*CVaR10, max(frequencyCounts)/1.9, 'CVaR10'); text(0.7*VaR10, max(frequencyCounts)/1.9, 'VaR10');
title('10-day Loss horizon')
xlabel('Loss')
ylabel('Number of Scenario')

%%%%%%%%%%%Part 2%%%%%%%%%%%
% Var for three stocks by using historical data
MSFT_loss = loss_1d(:,1);
AAPL_loss = loss_1d(:,9);
IBM_loss = loss_1d(:,10);

MSFT_loss = sort(MSFT_loss *100);
AAPL_loss = sort(AAPL_loss *200);
IBM_loss = sort(IBM_loss * 500);

VaR1_MSFT  = MSFT_loss(ceil(Ns_1 * alf));
VaR1_AAPL  = AAPL_loss(ceil(Ns_1 * alf));
VaR1_IBM  = IBM_loss(ceil(Ns_1 * alf));

SumVar = VaR1_MSFT + VaR1_AAPL + VaR1_IBM;
Diff_Var = VaR1 - SumVar;
fprintf('\npart 2\n')
fprintf('Historical 1-day MSFT VaR %4.1f%% = $%6.2f\n', 100*alf, VaR1_MSFT)
fprintf('Historical 1-day AAPL VaR %4.1f%% = $%6.2f\n', 100*alf, VaR1_AAPL)
fprintf('Historical 1-day IBM VaR %4.1f%% = $%6.2f\n', 100*alf, VaR1_IBM)
fprintf('Historical 1-day Simple Sum of VaR is $%6.2f% and the portfolio would give $%6.2f\n', SumVar, VaR1)
fprintf('Historical 1-day Difference in sum and Portfolio at %4.1f% is = $%6.2f\n', 100*alf, Diff_Var)

% Question 2
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio128\cplex\matlab\x64_win64');
% Annual risk-free rate for years 2015-2016 is 2.5%
r_rf = 0.025;

% Initial portfolio weights
init_positions = [5000 950 2000 0 0 0 0 2000 3000 1500 0 0 0 0 0 0 1001 0 0 0]';
init_value = data_prices(day_ind_end+1,:) * init_positions;
w_init = (data_prices(day_ind_end+1,:) .* init_positions')' / init_value;

% Max Sharpe Ratio portfolio weights
w_Sharpe = [ 0 0 0 0 0 0 0 0.385948690661642 0.172970428625544 0 0 0 0 0 0.003409676869715 0.260942060896445 0 0.185966939781285 0 0]';

% Equal Risk Contribution portfolio weights
w_ERC = [0.049946771209069 0.049951626261681 0.049955739901370 0.049998404150207 0.050000297368719 0.050004255546315 0.050006307026730 0.050007308995726 0.050010525832832 0.050013840015521 0.050014404492514 0.050015932843104 0.050016630302524 0.050017212457105 0.050017600497611 0.050017998351827 0.050018997074443 0.050019598350121 0.050019778113513 0.049946771209069]';

%%%%% Insert your code here 
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
w_eq = (1/n)*ones(n,1);
% Equally weighted portfolio return and variance
ret_equal_weight = dot(mu, w_eq);
var_equal_weight = w_eq' * Q * w_eq;

%Initial Portfolio 
ret_init_portf = dot(mu, w_init);
var_init_portf = w_init' * Q * w_init;

%Max Sharpe Ratio
ret_Sharpe = dot(mu, w_Sharpe);
var_Sharpe = w_Sharpe' * Q * w_Sharpe;

%Risk Free Asset
ret_rf = r_rf/252 ;
var_rf = 0;

%Eqully risk contribution
ret_ERC = dot(mu, w_ERC);
var_ERC = w_ERC' * Q * w_ERC;

%Levergeded equaly weight 
w_leqr = 2*w_ERC;
ret_LERC = 2 * ret_ERC - ret_rf;
var_LERC = w_leqr' * Q * w_leqr;

% Plot for Question 2, Part 1
% figure(3);
figure(3);
plot(sqrt(var_front), ret_front, 'k-', 'LineWidth', 3)
hold on;
plot(sqrt(var_minVar), ret_minVar, 'rd', 'MarkerSize', 6)
hold on;
plot(sqrt(var_maxRet), ret_maxRet, 'ms', 'MarkerSize', 6)
hold on;
plot(sqrt(var_equal_weight), ret_equal_weight, 'v', 'MarkerSize', 6)
hold on;
plot(sqrt(var_init_portf), ret_init_portf, 'h', 'MarkerSize', 6)
hold on;
plot(sqrt(var_Sharpe), ret_Sharpe, 'p', 'MarkerSize', 6)
hold on;
plot(sqrt(var_rf), ret_rf, 'bx', 'MarkerSize', 6)
hold on;
plot(sqrt(var_ERC), ret_ERC, 'm^', 'MarkerSize', 6)
hold on;
plot(sqrt(var_LERC), ret_LERC, 'ch', 'MarkerSize', 6)
hold on;
plot(sqrt(diag(Q)), mu, 'b.', 'MarkerSize', 18)
xlabel('Standard Deviation');
ylabel('Return');
title('Effient Frontier and Different Portfolio')
legend('Efficient Frontier', 'Minimum Variance Portfolio', 'Maximum Return Portfolio', 'Equal-weighted Portfolio', 'Initial portfolio', 'Maximum Sharpe Ratio portfolio', 'Risk-free Portfolio', 'Equal risk contribution Portfolio', 'Leverage equal risk contribution portfolio', 'Individual stocks', 'Location', 'SouthEast')
hold off

% Plot for Question 2, Part 2
% figure(4);
for j =1:1000
    rdm=rand(1,20);
    w_rand = rdm/sum(rdm);
    %calculate variance and return for each period 
    var_rd(j) = w_rand*Q*w_rand';
    ret_rd(j) = mu'* w_rand';
end
figure(4);
set(gcf, 'color', 'white');
plot(sqrt(var_front), ret_front, 'k-', 'LineWidth', 3)
hold on;
plot(sqrt(diag(Q)), mu, 'b.', 'MarkerSize', 18)
hold on;
for j=1:length(var_rd)
    plot(sqrt(var_rd(j)), ret_rd(j), 'r.', 'MarkerSize', 10)
end
hold on 
xlabel('Standard Deviation');
ylabel('Expected Return');
title('Efficient Frontier and Random Weighted Assets');
legend('Effiencient frontier','Individual Stock','Random Assets');
hold off
