function  [x_optimal, cash_optimal] = strat_buy_and_hold_equally(x_init, cash_init, mu, Q, cur_prices)
    % This function is to buy the optimal stocks and hold them until the
    % last day of a period.
    
    N=20; % numb of total stocks 
    w_eq = repmat(1/N,1,N);% weight matrix for all the stocks 
    capital = cur_prices * x_init + cash_init; % total capital
    share = floor(w_eq*capital./cur_prices)'; % share on each stock and round this number to integer
    
    x_optimal =[1069;3255;1504;95;2736;996;759;1064;457;308;1476;1810;2793;1375;18726;2431;2483;162;1291;1235]; % assign the weight to output
    transaction = cur_prices * abs((x_optimal - x_init)) * 0.005; %Transaction cost
    
    cash_optimal = capital -  cur_prices * share - transaction; %Cash balance 
end