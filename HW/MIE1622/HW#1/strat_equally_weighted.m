function  [x_optimal cash_optimal] = strat_equally_weighted(x_init, cash_init, mu, Q, cur_prices)
    %Equally weighted stragety distributes total asset into n portions into each stock
    %in our case n=20 
    
    N=20; % numb of total stocks 
    w_eq = repmat(1/N,1,N);% weight matrix for all the stocks 
    capital = cur_prices * x_init + cash_init; % total capital
    share = floor(w_eq*capital./cur_prices)'; % share on each stock and round this number to integer
    
    x_optimal = share; % assign the weight to output 
    
    transaction = cur_prices * abs((x_optimal - x_init)) * 0.005; %Transaction cost
    
    cash_optimal = capital -  cur_prices * share - transaction; %Cash balance 
    
    %check whether cash account has negative balance recalculation will be
    %done in the main function 
    if cash_optimal < 0
        fprintf('ALERT: "Equally Weighted Portfolio" strategy - CASH ACCOUNT BELOW ZERO');
    end
    
end
