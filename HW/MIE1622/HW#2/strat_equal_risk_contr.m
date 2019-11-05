function  [x_optimal cash_optimal] = strat_equal_risk_contr(x_init, cash_init, mu, var, cur_prices, rf, init_value)

    global A_ineq A_eq Q
    
    % number of stocks involved in previous period
    n = length(x_init);
    A_eq = ones(1,n);
    b_eq = 1; 
    Q = var;      % global variance of Q
    
    % Inequality constraints
    A_ineq = [];
    b_ineql = [];
    b_inequ = [];
           
    % Define initial portfolio 
    pot_value1 = cur_prices * x_init+cash_init;    % calcualte last period portfolio values
    w0 = cur_prices'.*x_init ./ pot_value1;
    %w0 = repmat(1.0/n, n, 1);      

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
    x_optimal = floor((pot_value1.*w_erc)./cur_prices');
    cash_optimal = pot_value1 - cur_prices*(x_optimal);
    
end