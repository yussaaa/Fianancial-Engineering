function  [x_optimal cash_optimal] = strat_robust_optim(x_init, cash_init, mu, var, cur_prices, r_rf, init_value)
    
    n = length(x_init);
    % Define initial portfolio 
    captial = cur_prices * x_init+cash_init;    % calcualte last period portfolio values
    w0 = cur_prices'.*x_init ./ captial;        % calcualte initial weightings
    
    global Q  
    Q = var;    % read variance value from main function 
    
    % Bounds on variables
    lb_rMV = zeros(n,1);
    ub_rMV = inf*ones(n,1);

    % Required portfolio robustness
    var_matr = diag(diag(Q));
    % Target portfolio return estimation error is return estimation error of 1/n portfolio
    rob_init = w0' * var_matr * w0; % return estimation error of initial portfolio
    rob_bnd = rob_init; % target return estimation error

    % Compute minimum variance portfolio
    cplex_minVar = Cplex('MinVar');
    cplex_minVar.DisplayFunc = []; % disable output to screen
    cplex_minVar.addCols(zeros(1,n)', [], lb_rMV, ub_rMV);
    cplex_minVar.addRows(1, ones(1,n), 1);
    cplex_minVar.Model.Q = 2*Q;
    cplex_minVar.Param.qpmethod.Cur = 6;
    cplex_minVar.solve();
    %cplex_minVar.Solution

    w_minVar = cplex_minVar.Solution.x; % asset weights
    ret_minVar = dot(mu, w_minVar); % Min Var portfolio return
    var_minVar = w_minVar' * Q * w_minVar; %  Min Var portfolio variance 
    rob_minVar = w_minVar' * var_matr * w_minVar;

    % Target portfolio return is return of minimum variance portfolio
    Portf_Retn = ret_minVar;

    % Objective function
    f_rMV  = zeros(n,1);
    % Constraints
    A_rMV  = sparse([  mu';...
                 ones(1,n)]);
    lhs_rMV = [Portf_Retn; 1];
    rhs_rMV = [inf; 1];
    % Initialize CPLEX environment
    cplex_rMV = Cplex('Robust_MV');
    cplex_rMV.DisplayFunc = []; % disable output to screen

    % Add objective function and variable bounds
    cplex_rMV.addCols(f_rMV, [], lb_rMV, ub_rMV);
    % Add constraints
    cplex_rMV.addRows(lhs_rMV, A_rMV, rhs_rMV);
    % Add quadratic objective
    cplex_rMV.Model.Q = 2*Q;
    % Add quadratic constraint on return estimation error (robustness constraint)
    Qq_rMV = var_matr;
    cplex_rMV.addQCs(zeros(size(f_rMV)), Qq_rMV, 'L', rob_bnd, {'qc_robust'});
    % Set CPLEX parameters
    cplex_rMV.Param.threads.Cur = 4;
    cplex_rMV.Param.timelimit.Cur = 60;
    cplex_rMV.Param.barrier.qcpconvergetol.Cur = 1e-12; % solution tolerance
    cplex_rMV.solve(); 
    %cplex_rMV.Solution 

    w_rMV = cplex_rMV.Solution.x;    
    
    % Round near-zero portfolio weights
    w_rMV_nonrnd = w_rMV;
    w_rMV(find(w_rMV<=1e-6)) = 0;
    w_rMV = w_rMV / sum(w_rMV);

    x_optimal = floor((captial.*w_rMV)./cur_prices');
    cash_optimal = captial - cur_prices*(x_optimal);
end