function  [x_optimal cash_optimal] = strat_lever_equal_risk_contr(x_init, cash_init, mu, var, cur_prices, r_rf, init_value1)

    global A_ineq A_eq
    global Q 
    global period;
%     global init_value2008;
%     global init_value;
    
    % number of stocks involved in previous period
    n = length(x_init);
    A_eq = ones(1,n);
    b_eq = 1; 
    Q = var;      % global variance of Q
    r_rff = r_rf/6;   % risk free rate per period
    
    % Inequality constraints
    A_ineq = [];
    b_ineql = [];
    b_inequ = [];
           
    % Define initial portfolio 
    capital = cur_prices * x_init + cash_init;    % calcualte last period portfolio values
    borrowed = capital;
    w0 = cur_prices'.*x_init ./ capital;
    %w0 = repmat(1.0/n, n, 1);      
    
    % Since this method borrows money in the first palce, we specify the
    % total asset and borrowed money
    if (period == 1)
        capital = (cur_prices * x_init+cash_init) + borrowed; 
    else 
        capital = (cur_prices * x_init+cash_init);   % hold/ not change it for the rest of periods 
    end

    options.lb = zeros(1,n);       % lower bounds on variables
    options.lu = ones (1,n);       % upper bounds on variables
    options.cl = [b_eq' b_ineql']; % lower bounds on constraints
    options.cu = [b_eq' b_inequ']; % upper bounds on constraints

    % Set the IPOPT options
    options.ipopt.jac_c_constant        = 'yes';
    options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.mu_strategy           = 'adaptive';
    options.ipopt.tol                   = 1e-10;
    options.ipopt.print_level = 0;

    % The callback functions
    funcs.objective         = @computeObjERC;
    funcs.constraints       = @computeConstraints;
    funcs.gradient          = @computeGradERC;
    funcs.jacobian          = @computeJacobian;
    funcs.jacobianstructure = @computeJacobian;

    % !!!! Function "computeGradERC" is just the placeholder
    % !!!! You need to compute the gradient yourself
    % Function "computeGradERC" revised. 

    % Run IPOPT
    [wsol info] = ipopt(w0',funcs,options);

    % Make solution a column vector
    if(size(wsol,1)==1)
        w_erc = wsol';
    else
        w_erc = wsol;
    end

    % Compute return, variance and risk contribution for the ERC portfolio
    ret_ERC = dot(mu, w_erc);
    var_ERC = w_erc'*Q*w_erc;
    RC_ERC = (w_erc .* ( Q*w_erc )) / sqrt(w_erc'*Q*w_erc);
    %disp(RC_ERC)
    
    % output
    x_optimal = floor((capital.*w_erc)./cur_prices');
    transaction = cur_prices * abs((x_optimal - x_init)) * 0.005;
    period_interest = borrowed*r_rff;
    cash_optimal = capital - cur_prices * x_optimal - transaction - period_interest;
    % Do not need to pay for the interest rate in first period
%     if (period == 1) 
%         cash_optimal = (capital - cur_prices*(share)); 
%     % Need to pay back the borrowed money in the first place
%     else
%         cash_optimal = (capital - cur_prices*(share))-init_value1*(r_rff); % final cash account with payback of interest
%     end 

%     x_optimal = share;
    
end
