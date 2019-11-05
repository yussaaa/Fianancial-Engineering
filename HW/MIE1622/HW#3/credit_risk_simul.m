clear all;
clc
format long;

Nout  = 100000; % number of out-of-sample scenarios
Nin   = 5000;   % number of in-sample scenarios
Ns    = 5;      % number of idiosyncratic scenarios for each systemic

C = 8;          % number of credit states

% Filename to save out-of-sample scenarios
filename_save_out  = 'scen_out';

% Read and parse instrument data
instr_data = dlmread('instrum_data.csv', ',');
instr_id   = instr_data(:,1);           % ID
driver     = instr_data(:,2);           % credit driver
beta       = instr_data(:,3);           % beta (sensitivity to credit driver)
recov_rate = instr_data(:,4);           % expected recovery rate
value      = instr_data(:,5);           % value
prob       = instr_data(:,6:6+C-1);     % credit-state migration probabilities (default to A)
exposure   = instr_data(:,6+C:6+2*C-1); % credit-state migration exposures (default to A)
retn       = instr_data(:,6+2*C);       % market returns

K = size(instr_data, 1); % number of  counterparties

% Read matrix of correlations for credit drivers
rho = dlmread('credit_driver_corr.csv', '\t');
sqrt_rho = (chol(rho))'; % Cholesky decomp of rho (for generating correlated Normal random numbers)

disp('======= Credit Risk Model with Credit-State Migrations =======')
disp('============== Monte Carlo Scenario Generation ===============')
disp(' ')
disp(' ')
disp([' Number of out-of-sample Monte Carlo scenarios = ' int2str(Nout)])
disp([' Number of in-sample Monte Carlo scenarios = ' int2str(Nin)])
disp([' Number of counterparties = ' int2str(K)])
disp(' ')

% Find credit-state for each counterparty
% 8 = AAA, 7 = AA, 6 = A, 5 = BBB, 4 = BB, 3 = B, 2 = CCC, 1 = default
[Ltemp, CS] = max(prob, [], 2); % Find the indeces of the maximum values of A 
clear Ltemp

% Account for default recoveries
exposure(:, 1) = (1-recov_rate) .* exposure(:, 1);

% Compute credit-state boundaries
CS_Bdry = norminv( cumsum(prob(:,1:C-1), 2) );

% -------- Insert your code here -------- %

if(~exist('scenarios_out.mat','file'))
    Losses_out = zeros(Nout,K);
    
    % -------- Insert your code here -------- %

    for s = 1:Nout
        % -------- Insert your code here -------- %
        x = normrnd(0,1,[1,50]);
        y = x*sqrt_rho; % Correlacted random number
        z = normrnd(0,1,[100,1]);
        loss = zeros(1,100);
        w=zeros(K,1);
        for t = 1:K
            w(t)= beta(t)*y(driver(t))+sqrt(1-beta(t).^2)'*z(t);
            for i = 1:6     % check credit boundary and assign exposures to each credit rating
                if (w(t)>= CS_Bdry(t,i)) && (w(t)<CS_Bdry(t,i+1))
                    loss(t)= exposure(t,i+1);         % calculate portfolio loss
                elseif (w(t)< CS_Bdry(t,1))
                    loss(t)= exposure(t,1);
                else
                    loss(t)=exposure(t,8);
                end
            end
        end
        Losses_out(s,:)=loss;       % Calculated out-of-sample losses (100000 x 100)
    end

    % Calculated out-of-sample losses (100000 x 100)
    % Losses_out

    save('scenarios_out', 'Losses_out')
else
    load('scenarios_out', 'Losses_out')
end

%% Normal approximation computed from out-of-sample scenarios
mu_l = mean(Losses_out)';
var_l = cov(Losses_out);

% Compute portfolio weights
portf_v = sum(value);     % portfolio value
w0{1} = value / portf_v;  % asset weights (portfolio 1)
w0{2} = ones(K, 1) / K;   % asset weights (portfolio 2)
x0{1} = (portf_v ./ value) .* w0{1};  % asset units (portfolio 1)
x0{2} = (portf_v ./ value) .* w0{2};  % asset units (portfolio 2)

% Quantile levels (99%, 99.9%)
alphas = [0.99 0.999];

Losses_out_wtd{1}=sort(Losses_out*x0{1}); % total loss for 10000 scenarios using port 1
Losses_out_wtd{2}=sort(Losses_out*x0{2}); % total loss for 10000 scenarios using port 2

% Compute VaR and CVaR (non-Normal and Normal) for 100000 scenarios
for(portN = 1:2)
    for(q=1:length(alphas))
        alf = alphas(q);
%         x = x0{portN};
%         w= w0{portN};
%         port_loss_out{portN} = sort(Losses_out * x0{portN});
%         portloss_mu = mean(port_loss);
%         port_loss_std = std(port_loss);

%         VaRout(portN,q)  = port_loss(ceil(Nout*alf));
%         VaRinN(portN,q)  = portloss_mu + norminv(alf,0,1)*port_loss_std;
%         CVaRout(portN,q) = (1/(Nout*(1-alf))) * ( (ceil(Nout*alf)-Nout*alf) * VaRout(portN,q) + sum(port_loss(ceil(Nout*alf)+1:Nout)) );
%         CVaRinN(portN,q) = portloss_mu + (normpdf(norminv(alf,0,1))/(1-alf))*port_loss_std;
%         % -------- Insert your code here -------- %        

        VaRout(portN,q)  = Losses_out_wtd{portN}(ceil(Nout*alf));
        VaRinN(portN,q)  = mean(Losses_out_wtd{portN}) + norminv(alf,0,1)*std(Losses_out_wtd{portN});
        CVaRout(portN,q) = (1/(Nout*(1-alf))) * ( (ceil(Nout*alf)-Nout*alf) * VaRout(portN,q)+ sum(Losses_out_wtd{portN}(ceil(Nout*alf)+1:Nout)));
        CVaRinN(portN,q) = mean(Losses_out_wtd{portN}) + (normpdf(norminv(alf,0,1))/(1-alf))*std(Losses_out_wtd{portN});
        % -------- Insert your code here -------- %  
    end
end

      


% Perform 100 trials
N_trials = 100;

for(tr=1:N_trials)
%% Monte Carlo approximation 1
    % -------- Insert your code here -------- %
    for s = 1:ceil(Nin/Ns) % systemic scenarios
        % -------- Insert your code here -------- %
        x1 = normrnd(0,1,[50,1]);
        y1 = sqrt_rho*x1;
%         y1(s,:) = (sqrt_rho*x1)';
%         z1 = normrnd(0,1,5,1); %%%method 2
        for si = 1:Ns % idiosyncratic scenarios for each systemic
            % -------- Insert your code here -------- %
            for t = 1:K
            z1 = normrnd(0,1);
            w1(t)= beta(t)*y1(driver(t))+sqrt(1-beta(t).^2)'*z1;
            ind_m1 = find(sort([CS_Bdry(t,:),w1(t)]) == w1(t));
            Losses_inMC1(Ns*(s-1)+si,t) = exposure(t,ind_m1);
            
% % % %             for i = 1:6     % check credit boundary and assign exposures to each credit rating
% % % %                 if (w(t)>= CS_Bdry(t,i)) && (w(t)<CS_Bdry(t,i+1))
% % % %                     loss(t)= exposure(t,i+1);         % calculate portfolio loss
% % % %                 elseif (w(t)< CS_Bdry(t,1))
% % % %                     loss(t)= exposure(t,1);
% % % %                 else
% % % %                     loss(t)=exposure(t,8);
% % % %                 end
% % % %             end
%                  Losses_mc1{tr}((s-1)*Ns+si,t)=loss1(t);       % Calculated out-of-sample losses (100000 x 100)
             
            end    
        end
    end
    
    % Calculated losses for MC1 approximation (5000 x 100)
    % Losses_inMC1
    
    %% Monte Carlo approximation 2
    
    % -------- Insert your code here -------- %
    
    for s = 1:Nin % systemic scenarios (1 idiosyncratic scenario for each systemic)
        % -------- Insert your code here -------- %
        x2 = normrnd(0,1,[1,50]);
        y2 = x2*sqrt_rho;
        z2 = normrnd(0,1,[1,100]);
        for t = 1:K
            w2(t)= beta(t)*y2(driver(t))+sqrt(1-beta(t).^2)*z2(t);
% % % % % % %             ind_m2 = find(sort([CS_Bdry(t,:),w2(t)]) == w2(t));
% % % % % % % %             Losses_inMC2{tr}(Ns*(s-1)+si,t) = exposure(t,ind_m2);
% % % % % % % %             %%%%%?????????????????????????????????????????????????
% % % % % % %             Losses_inMC2(s,t) = exposure(t,ind_m2);
            for i = 1:6     % check credit boundary and assign exposures to each credit rating
                if (w2(t)>= CS_Bdry(t,i)) && (w2(t)<CS_Bdry(t,i+1))
                    loss2(t)= exposure(t,i+1);         % calculate portfolio loss
                elseif (w2(t)< CS_Bdry(t,1))
                    loss2(t)= exposure(t,1);
                else
                    loss2(t)=exposure(t,8);
                end
            end
        end
        Losses_inMC2(s,:)=loss2;       % Calculated losses (5000 x 100)
    
    end  
    

    % Calculated losses for MC2 approximation (5000 x 100)
    % Losses_inMC2
    
    %% Compute VaR and CVaR
    for(portN = 1:2)
        for(q=1:length(alphas))
            alf = alphas(q);
            % -------- Insert your code here -------- %       
            
%             x = x0{portN};
            w = w0{portN};
            % Compute portfolio loss 
            portf_loss_inMC1{portN} = sort(Losses_inMC1 * x0{portN});
            portf_loss_inMC2{portN} = sort(Losses_inMC2 * x0{portN});
            mu_MCl = mean(Losses_inMC1)';
            std_MCl = std(Losses_inMC1)';
            var_MCl = cov(Losses_inMC1);
            mu_MC2 = mean(Losses_inMC2)';
            std_MC2 = std(Losses_inMC2)';
            var_MC2 = cov(Losses_inMC2);
            % Compute portfolio mean loss mu_p_MC1 and portfolio standard deviation of losses sigma_p_MC1
            mu_p_MC1{portN} = mu_MCl'*x0{portN};
            sigma_p_MC1{portN} = std(portf_loss_inMC1{portN});
            % Compute portfolio mean loss mu_p_MC2 and portfolio standard deviation of losses sigma_p_MC2
            mu_p_MC2{portN} = mu_MC2'*x0{portN};
            sigma_p_MC2{portN} = std(portf_loss_inMC2{portN});
            % Compute VaR and CVaR for the current trial
            VaRinMC1{portN,q}(tr) = portf_loss_inMC1{portN}(ceil(Nin*alf));
            VaRinMC2{portN,q}(tr) = portf_loss_inMC2{portN}(ceil(Nin*alf));
            VaRinN1{portN,q}(tr) = mu_p_MC1{portN} + norminv(alf,0,1)*sigma_p_MC1{portN};
            VaRinN2{portN,q}(tr) = mu_p_MC2{portN} + norminv(alf,0,1)*sigma_p_MC2{portN};
            CVaRinMC1{portN,q}(tr) = (1/(Nin*(1-alf))) * ( (ceil(Nin*alf)-Nin*alf) * VaRinMC1{portN,q}(tr) + sum(portf_loss_inMC1{portN}(ceil(Nin*alf)+1:Nin)) );
            CVaRinMC2{portN,q}(tr) = (1/(Nin*(1-alf))) * ( (ceil(Nin*alf)-Nin*alf) * VaRinMC2{portN,q}(tr) + sum(portf_loss_inMC2{portN}(ceil(Nin*alf)+1:Nin)) );
            CVaRinN1{portN,q}(tr) = mu_p_MC1{portN} + (normpdf(norminv(alf,0,1))/(1-alf))*sigma_p_MC1{portN};
            CVaRinN2{portN,q}(tr) = mu_p_MC2{portN} + (normpdf(norminv(alf,0,1))/(1-alf))*sigma_p_MC2{portN};
        
            
            
%             Compute portfolio loss 
%             portf_loss_inMC1 = sort(Losses_inMC1 * x0{portN});
%             portf_loss_inMC2{tr,portN} = sort(Losses_inMC2 * x0{portN});
%             mu_MCl = mean(Losses_inMC1)';
%             var_MCl = cov(Losses_inMC1);
%             mu_MC2 = mean(Losses_inMC2)';
%             var_MC2 = cov(Losses_inMC2);
% %             Compute portfolio mean loss mu_p_MC1 and portfolio standard deviation of losses sigma_p_MC1
%             mu_p_MC1=mu_MCl'*x0{portN};
%             sigma_p_MC1=std(portf_loss_inMC1{tr,portN});
% %             Compute portfolio mean loss mu_p_MC2 and portfolio standard deviation of losses sigma_p_MC2
%             mu_p_MC2=mu_MC2'*x0{portN};
%             sigma_p_MC2=std(portf_loss_inMC2{tr,portN});
%        
% % %%%%%%%%%%%%%Compute VaR and CVaR for the current trial
%             VaRinMC1{portN,q}(tr) = portf_loss_inMC1{tr,portN}(ceil(Nin*alf));
%             VaRinMC2{portN,q}(tr) = portf_loss_inMC2{tr,portN}(ceil(Nin*alf));
%             VaRinN1{portN,q}(tr) = mu_p_MC1 + norminv(alf,0,1)*sigma_p_MC1;
%             VaRinN2{portN,q}(tr) = mu_p_MC2 + norminv(alf,0,1)*sigma_p_MC2;
%             CVaRinMC1{portN,q}(tr) = (1/(Nin*(1-alf))) * ( (ceil(Nin*alf)-Nin*alf) * VaRinMC1{portN,q}(tr) + sum(portf_loss_inMC1(ceil(Nin*alf)+1:Nin)) );
%             CVaRinMC2{portN,q}(tr) = (1/(Nin*(1-alf))) * ( (ceil(Nin*alf)-Nin*alf) * VaRinMC2{portN,q}(tr) + sum(portf_loss_inMC2(ceil(Nin*alf)+1:Nin)) );
%             CVaRinN1{portN,q}(tr) = mu_p_MC1 + (normpdf(norminv(alf,0,1))/(1-alf))*std_p_MC1;
%             CVaRinN2{portN,q}(tr) = mu_p_MC2 + (normpdf(norminv(alf,0,1))/(1-alf))*std_p_MC2;
%             CVaRinN2{portN,q}(tr) = ...
            % -------- Insert your code here -------- %
        end
    end
end

% Display portfolio VaR and CVaR
for(portN = 1:2)
fprintf('\nPortfolio %d:\n\n', portN)    
 for(q=1:length(alphas))
    alf = alphas(q);
    fprintf('Out-of-sample: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, VaRout(portN,q), 100*alf, CVaRout(portN,q))
    fprintf('In-sample MC1: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, mean(VaRinMC1{portN,q}), 100*alf, mean(CVaRinMC1{portN,q}))
    fprintf('In-sample MC2: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, mean(VaRinMC2{portN,q}), 100*alf, mean(CVaRinMC2{portN,q}))
    fprintf(' In-sample No: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, VaRinN(portN,q), 100*alf, CVaRinN(portN,q))
    fprintf(' In-sample N1: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, mean(VaRinN1{portN,q}), 100*alf, mean(CVaRinN1{portN,q}))
    fprintf(' In-sample N2: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n\n', 100*alf, mean(VaRinN2{portN,q}), 100*alf, mean(CVaRinN2{portN,q}))
 end
end

%% Plot results
% figure(1);
% -------- Insert your code here -------- %
figure(1);
set(gcf, 'color', 'white');
[frequencyCounts, binLocations] = hist(Losses_out_wtd{1}, 100);
bar(binLocations, frequencyCounts);
hold on;
line([VaRout(1,1) VaRout(1,1)], [0 max(frequencyCounts)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
hold on;
line([VaRout(1,2) VaRout(1,2)], [0 max(frequencyCounts)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
hold on;
line([CVaRout(1,1) CVaRout(1,1)], [0 max(frequencyCounts)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')
hold on;
line([CVaRout(1,2) CVaRout(1,2)], [0 max(frequencyCounts)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% line([VaRinN(1,1) VaRinN(1,1)], [0 max(frequencyCounts)], 'Color', 'y', 'LineWidth', 0.5, 'LineStyle', '--')
% line([VaRinN(1,2) VaRinN(1,2)], [0 max(frequencyCounts)], 'Color', 'y', 'LineWidth', 0.5, 'LineStyle', '--')
% line([CVaRinN(1,1) CVaRinN(1,1)], [0 max(frequencyCounts)], 'Color', 'g', 'LineWidth', 0.5, 'LineStyle', '--')
% line([CVaRinN(1,2) CVaRinN(1,2)], [0 max(frequencyCounts)], 'Color', 'g', 'LineWidth', 0.5, 'LineStyle', '--')

normf = ( 1/(std(Losses_out_wtd{1})*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(Losses_out_wtd{1}))/std(Losses_out_wtd{1})).^2 );
normf = normf * sum(frequencyCounts)/sum(normf);
plot(binLocations, normf, 'r', 'LineWidth', 1); hold on;

text(1.2*VaRout(1,1), max(frequencyCounts)/1.5, 'VaR99%')
text(0.7*CVaRout(1,1), max(frequencyCounts)/1.3, 'CVaR99%')
text(1.1*VaRout(1,2), max(frequencyCounts)/1.4, 'VaR99.9%')
text(0.7*CVaRout(1,2), max(frequencyCounts)/1.1, 'CVaR99.9%')

% text(1.2*VaRinN(1,1), max(frequencyCounts)/1.5, 'VaRn99%')
% text(0.7*CVaRinN(1,1), max(frequencyCounts)/1.3, 'CVaRn99%')
% text(1.1*VaRinN(1,2), max(frequencyCounts)/1.4, 'VaRn99.9%')
% text(0.7*CVaRinN(1,2), max(frequencyCounts)/1.1, 'CVaRn99.9%')

hold off;
title('Out-of-sample loss distribution, portfolio 1 - VaR and CVaR values')
xlabel('Loss $')
ylabel('Frequency')
% figure(2);
% -------- Insert your code here -------- 

% figure 2: out-of-sample distribution, portfolio 2 
figure(2);
set(gcf, 'color', 'white');
[frequencyCounts1, binLocations1] = hist(Losses_out_wtd{2}, 100);
bar(binLocations1, frequencyCounts1);
hold on;
line([VaRout(2,1) VaRout(2,1)], [0 max(frequencyCounts1)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
% hold on;
line([VaRout(2,2) VaRout(2,2)], [0 max(frequencyCounts1)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
% hold on;
line([CVaRout(2,1) CVaRout(2,1)], [0 max(frequencyCounts1)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')
% hold on;
line([CVaRout(2,2) CVaRout(2,2)], [0 max(frequencyCounts1)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')
% hold on;
normf = ( 1/(std(Losses_out_wtd{2})*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(Losses_out_wtd{2}))/std(Losses_out_wtd{2})).^2 );
normf = normf * sum(frequencyCounts)/sum(normf);
plot(binLocations, normf, 'r', 'LineWidth', 1); %hold on;

text(1.2*VaRout(2,1), max(frequencyCounts1)/1.5, 'VaR99%')
text(0.7*CVaRout(2,1), max(frequencyCounts1)/1.3, 'CVaR99%')
text(1.1*VaRout(2,2), max(frequencyCounts1)/1.4, 'VaR99.9%')
text(0.7*CVaRout(2,2), max(frequencyCounts1)/1.1, 'CVaR99.9%')
title('Out-of-sample loss distribution, portfolio 2 - VaR and CVaR values')
xlabel('Loss $')
ylabel('Frequency')
hold off;


%% MC1
% figure 3: MC 1, portfolio 1
figure(3);
set(gcf, 'color', 'white');
[frequencyCounts2, binLocations2] = hist(portf_loss_inMC1{1}, 100);
bar(binLocations2, frequencyCounts2);
hold on;
line([mean(VaRinMC1{1,1}) mean(VaRinMC1{1,1})], [0 max(frequencyCounts2)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
hold on;
line([mean(VaRinMC1{1,2}) mean(VaRinMC1{1,2})], [0 max(frequencyCounts2)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
hold on;
line([mean(CVaRinMC1{1,1}) mean(CVaRinMC1{1,1})], [0 max(frequencyCounts2)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')
hold on;
line([mean(CVaRinMC1{1,2}) mean(CVaRinMC1{1,2})], [0 max(frequencyCounts2)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')
hold on;
text(0.98*mean(VaRinMC1{1,1}), max(frequencyCounts2)/1.2, 'VaR99%')
text(0.98*mean(CVaRinMC1{1,1}), max(frequencyCounts2)/1.2, 'CVaR99%')
text(0.98*mean(VaRinMC1{1,2}), max(frequencyCounts2)/1.5, 'VaR99.9%')
text(0.98*mean(CVaRinMC1{1,2}), max(frequencyCounts2)/1.5, 'CVaR99.9%')

normf = ( 1/(std(portf_loss_inMC1{1})*sqrt(2*pi)) ) * exp( -0.5*((binLocations2-mean(portf_loss_inMC1{1}))/std(portf_loss_inMC1{1})).^2 );
normf = normf * sum(frequencyCounts2)/sum(normf);
plot(binLocations2, normf, 'r', 'LineWidth', 1);

title('MC1 loss distribution, portfolio 1 - VaR and CVaR values')
xlabel('Loss $')
ylabel('Frequency')
hold off;


clear frequencyCounts
clear binLocations
clear normf
% figure 4: MC 1, portfolio 2
figure(4);
set(gcf, 'color', 'white');
[frequencyCounts, binLocations] = hist(portf_loss_inMC1{2}, 100);
bar(binLocations, frequencyCounts);
hold on;
line([mean(VaRinMC1{2,1}) mean(VaRinMC1{2,1})], [0 max(frequencyCounts)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
hold on;
line([mean(VaRinMC1{2,2}) mean(VaRinMC1{2,2})], [0 max(frequencyCounts)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
hold on;
line([mean(CVaRinMC1{2,1}) mean(CVaRinMC1{2,1})], [0 max(frequencyCounts)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')
hold on;
line([mean(CVaRinMC1{2,2}) mean(CVaRinMC1{2,2})], [0 max(frequencyCounts)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')

text(0.9*mean(VaRinMC1{2,1}), max(frequencyCounts)/1.1, 'VaR99%')
text(0.9*mean(CVaRinMC1{2,1}), max(frequencyCounts)/1.1, 'CVaR99%')
text(0.9*mean(VaRinMC1{2,2}), max(frequencyCounts)/1.4, 'VaR99.9%')
text(0.9*mean(CVaRinMC1{2,2}), max(frequencyCounts)/1.4, 'CVaR99.9%')

title('MC1 loss distribution, portfolio 2 - VaR and CVaR values')
xlabel('Loss $')
ylabel('Frequency')

normf = ( 1/(std(portf_loss_inMC1{2})*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(portf_loss_inMC1{2}))/std(portf_loss_inMC1{2})).^2 );
normf = normf * sum(frequencyCounts)/sum(normf);
plot(binLocations, normf, 'r', 'LineWidth', 1); %hold on;
hold off;

clear frequencyCounts
clear binLocations
clear normf
%% MC2
% figure 5: MC 2, portfolio 1
figure(5);
set(gcf, 'color', 'white');
[frequencyCounts, binLocations] = hist(portf_loss_inMC2{1}, 100);
bar(binLocations, frequencyCounts);
hold on;
line([mean(VaRinMC2{1,1}) mean(VaRinMC2{1,1})], [0 max(frequencyCounts)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
hold on;
line([mean(VaRinMC2{1,2}) mean(VaRinMC2{1,2})], [0 max(frequencyCounts)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
hold on;
line([mean(CVaRinMC2{1,1}) mean(CVaRinMC2{1,1})], [0 max(frequencyCounts)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')
hold on;
line([mean(CVaRinMC2{1,2}) mean(CVaRinMC2{1,2})], [0 max(frequencyCounts)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')

text(1.2*mean(VaRinMC2{1,1}), max(frequencyCounts)/1.5, 'VaR99%')
text(0.7*mean(CVaRinMC2{1,1}), max(frequencyCounts)/1.3, 'CVaR99%')
text(1.1*mean(VaRinMC2{1,2}), max(frequencyCounts)/1.4, 'VaR99.9%')
text(0.7*mean(CVaRinMC2{1,2}), max(frequencyCounts)/1.1, 'CVaR99.9%')
title('MC2 loss distribution, portfolio 1 - VaR and CVaR values')
xlabel('Loss $')
ylabel('Frequency')
normf = ( 1/(std(portf_loss_inMC2{1})*sqrt(2*pi)) ) * exp( -0.5*((binLocations2-mean(portf_loss_inMC2{1}))/std(portf_loss_inMC2{1})).^2 );
normf = normf * sum(frequencyCounts2)/sum(normf);
plot(binLocations2, normf, 'r', 'LineWidth', 1);
hold off;

clear frequencyCounts
clear binLocations
clear normf
% figure 6: MC 2, portfolio 2
figure(6);
set(gcf, 'color', 'white');
[frequencyCounts, binLocations] = hist(portf_loss_inMC2{2}, 100);
bar(binLocations, frequencyCounts);
hold on;
line([mean(VaRinMC2{2,1}) mean(VaRinMC2{2,1})], [0 max(frequencyCounts)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
hold on;
line([mean(VaRinMC2{2,2}) mean(VaRinMC2{2,2})], [0 max(frequencyCounts)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
hold on;
line([mean(CVaRinMC2{2,1}) mean(CVaRinMC2{2,1})], [0 max(frequencyCounts)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')
hold on;
line([mean(CVaRinMC2{2,2}) mean(CVaRinMC2{2,2})], [0 max(frequencyCounts)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')

text(1.2*mean(VaRinMC2{2,1}), max(frequencyCounts)/1.5, 'VaR99%')
text(0.7*mean(CVaRinMC2{2,1}), max(frequencyCounts)/1.3, 'CVaR99%')
text(1.1*mean(VaRinMC2{2,2}), max(frequencyCounts)/1.4, 'VaR99.9%')
text(0.7*mean(CVaRinMC2{2,2}), max(frequencyCounts)/1.1, 'CVaR99.9%')
title('MC2 loss distribution, portfolio 2 - VaR and CVaR values')
xlabel('Loss $')
ylabel('Frequency')

normf = ( 1/(std(portf_loss_inMC2{2})*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(portf_loss_inMC2{2}))/std(portf_loss_inMC2{2})).^2 );
normf = normf * sum(frequencyCounts)/sum(normf);
plot(binLocations, normf, 'r', 'LineWidth', 1); %hold on;
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%% Normal %%%%%%%%%%%%%%%%%%%%%%%
%% In-sample No, portfolio 1
%%% figure 7 
clear frequencyCounts
clear binLocations
clear normf
figure(7);
set(gcf, 'color', 'white');
[frequencyCounts, binLocations] = hist(Losses_out_wtd{1}, 100);
bar(binLocations, frequencyCounts);
hold on;
% plot normal distribution curve on
normf = ( 1/(std(Losses_out_wtd{1})*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(Losses_out_wtd{1}))/std(Losses_out_wtd{1})).^2 );
normf = normf * sum(frequencyCounts)/sum(normf);
plot(binLocations, normf, 'r', 'LineWidth', 2);
hold on;
line([VaRinN(1,1) VaRinN(1,1)], [0 max(frequencyCounts)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
line([VaRinN(1,2) VaRinN(1,2)], [0 max(frequencyCounts)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
line([CVaRinN(1,1) CVaRinN(1,1)], [0 max(frequencyCounts)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')
line([CVaRinN(1,2) CVaRinN(1,2)], [0 max(frequencyCounts)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')
hold off;
text(1.2*VaRinN(1,1), max(frequencyCounts)/1.5, 'VaR99%')
text(0.7*CVaRinN(1,1), max(frequencyCounts)/1.3, 'CVaR99%')
text(1.1*VaRinN(1,2), max(frequencyCounts)/1.4, 'VaR99.9%')
text(0.7*CVaRinN(1,2), max(frequencyCounts)/1.1, 'CVaR99.9%')
title('Out-of-sample loss distribution(Normal), portfolio 1 - VaR and CVaR values')
xlabel('Loss $')
ylabel('Frequency')

%%% In-sample No, portfolio 2
clear frequencyCounts
clear binLocations
clear normf
figure(8);
set(gcf, 'color', 'white');
[frequencyCounts, binLocations] = hist(Losses_out_wtd{2}, 100);
bar(binLocations, frequencyCounts);
hold on;
% plot normal distribution curve on
normf = ( 1/(std(Losses_out_wtd{2})*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(Losses_out_wtd{2}))/std(Losses_out_wtd{2})).^2 );
normf = normf * sum(frequencyCounts)/sum(normf);
plot(binLocations, normf, 'r', 'LineWidth', 2);
hold on;
line([VaRinN(2,1) VaRinN(2,1)], [0 max(frequencyCounts)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
hold on;
line([VaRinN(2,2) VaRinN(2,2)], [0 max(frequencyCounts)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
hold on;
line([CVaRinN(2,1) CVaRinN(2,1)], [0 max(frequencyCounts)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')
hold on;
line([CVaRinN(2,2) CVaRinN(2,2)], [0 max(frequencyCounts)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')
hold off;
text(1.2*VaRinN(2,1), max(frequencyCounts)/1.5, 'VaR99%')
text(0.7*CVaRinN(2,1), max(frequencyCounts)/1.3, 'CVaR99%')
text(1.1*VaRinN(2,2), max(frequencyCounts)/1.4, 'VaR99.9%')
text(0.7*CVaRinN(2,2), max(frequencyCounts)/1.1, 'CVaR99.9%')
title('Out-of-sample loss distribution(Normal), portfolio 2 - VaR and CVaR values')
xlabel('Loss $')
ylabel('Frequency')


clear frequencyCounts
clear binLocations
clear normf

%% Normal

% figure 9: MC 1, portfolio 1, normal: mean(VaRinN1{portN,q})
figure(9);
set(gcf, 'color', 'white');
[frequencyCounts, binLocations] = hist(portf_loss_inMC1{1}, 100);
bar(binLocations, frequencyCounts);
hold on;
% plot normal distribution curve on
normf = ( 1/(std(portf_loss_inMC1{1})*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(portf_loss_inMC1{1}))/std(portf_loss_inMC1{1})).^2 );
normf = normf * sum(frequencyCounts)/sum(normf);
plot(binLocations, normf, 'r', 'LineWidth', 2);
hold on;
line([mean(VaRinN1{1,1}) mean(VaRinN1{1,1})], [0 max(frequencyCounts)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
hold on;
line([mean(VaRinN1{1,2}) mean(VaRinN1{1,2})], [0 max(frequencyCounts)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
hold on;
line([mean(CVaRinN1{1,1}) mean(CVaRinN1{1,1})], [0 max(frequencyCounts)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')
hold on;
line([mean(CVaRinN1{1,2}) mean(CVaRinN1{1,2})], [0 max(frequencyCounts)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')
hold off;
text(1.2*mean(VaRinN1{2,1}), max(frequencyCounts)/1.5, 'VaR99%')
text(0.7*mean(CVaRinN1{2,1}), max(frequencyCounts)/1.3, 'CVaR99%')
text(1.1*mean(VaRinN1{2,2}), max(frequencyCounts)/1.4, 'VaR99.9%')
text(0.7*mean(CVaRinN1{2,2}), max(frequencyCounts)/1.1, 'CVaR99.9%')
title('MC1 loss distribution, portfolio 1(normal) - VaR and CVaR values')
xlabel('Loss $')
ylabel('Frequency')


clear frequencyCounts
clear binLocations
clear normf
% figure 10: MC 1, portfolio 2, normal: mean(VaRinN1{portN,q})
figure(10);
set(gcf, 'color', 'white');
[frequencyCounts, binLocations] = hist(portf_loss_inMC1{2}, 100);
bar(binLocations, frequencyCounts);
hold on;
% plot normal distribution curve on
normf = ( 1/(std(portf_loss_inMC1{2})*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(portf_loss_inMC1{2}))/std(portf_loss_inMC1{2})).^2 );
normf = normf * sum(frequencyCounts)/sum(normf);
plot(binLocations, normf, 'r', 'LineWidth', 2);
hold on;
line([mean(VaRinN1{2,1}) mean(VaRinN1{2,1})], [0 max(frequencyCounts)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
hold on;
line([mean(VaRinN1{2,2}) mean(VaRinN1{2,2})], [0 max(frequencyCounts)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
hold on;
line([mean(CVaRinN1{2,1}) mean(CVaRinN1{2,1})], [0 max(frequencyCounts)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')
hold on;
line([mean(CVaRinN1{2,2}) mean(CVaRinN1{2,2})], [0 max(frequencyCounts)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')
hold off;
text(1.2*mean(VaRinN1{2,1}), max(frequencyCounts)/1.5, 'VaR99%')
text(0.7*mean(CVaRinN1{2,1}), max(frequencyCounts)/1.3, 'CVaR99%')
text(1.1*mean(VaRinN1{2,2}), max(frequencyCounts)/1.4, 'VaR99.9%')
text(0.7*mean(CVaRinN1{2,2}), max(frequencyCounts)/1.1, 'CVaR99.9%')
title('MC1 loss distribution, portfolio 2(normal) - VaR and CVaR values')
xlabel('Loss $')
ylabel('Frequency')


clear frequencyCounts
clear binLocations
clear normf
%% figure 11: MC 2, portfolio 1, normal
figure(11);
set(gcf, 'color', 'white');
[frequencyCounts, binLocations] = hist(portf_loss_inMC2{1}, 100);
bar(binLocations, frequencyCounts);
hold on;
% plot normal distribution curve on
normf = ( 1/(std(portf_loss_inMC2{1})*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(portf_loss_inMC2{1}))/std(portf_loss_inMC2{1})).^2 );
normf = normf * sum(frequencyCounts)/sum(normf);
plot(binLocations, normf, 'r', 'LineWidth', 2);
hold on;
line([mean(VaRinN2{1,1}) mean(VaRinN2{1,1})], [0 max(frequencyCounts)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
hold on;
line([mean(VaRinN2{1,2}) mean(VaRinN2{1,2})], [0 max(frequencyCounts)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
hold on;
line([mean(CVaRinN2{1,1}) mean(CVaRinN2{1,1})], [0 max(frequencyCounts)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')
hold on;
line([mean(CVaRinN2{1,2}) mean(CVaRinN2{1,2})], [0 max(frequencyCounts)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')
hold off;
text(1.2*mean(VaRinN2{1,1}), max(frequencyCounts)/1.5, 'VaR99%')
text(0.7*mean(CVaRinN2{1,1}), max(frequencyCounts)/1.3, 'CVaR99%')
text(1.1*mean(VaRinN2{1,2}), max(frequencyCounts)/1.4, 'VaR99.9%')
text(0.7*mean(CVaRinN2{1,2}), max(frequencyCounts)/1.1, 'CVaR99.9%')
title('MC2 loss distribution, portfolio 1(normal) - VaR and CVaR values')
xlabel('Loss $')
ylabel('Frequency')


clear frequencyCounts
clear binLocations
clear normf
% figure 12: MC 2, portfolio 2, normal
figure(12);
set(gcf, 'color', 'white');
[frequencyCounts, binLocations] = hist(portf_loss_inMC2{2}, 100);
bar(binLocations, frequencyCounts);
hold on;
% plot normal distribution curve on
normf = ( 1/(std(portf_loss_inMC2{2})*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(portf_loss_inMC2{2}))/std(portf_loss_inMC2{2})).^2 );
normf = normf * sum(frequencyCounts)/sum(normf);
plot(binLocations, normf, 'r', 'LineWidth', 2);
hold on;
line([mean(VaRinN2{2,1}) mean(VaRinN2{2,1})], [0 max(frequencyCounts)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
hold on;
line([mean(VaRinN2{2,2}) mean(VaRinN2{2,2})], [0 max(frequencyCounts)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
hold on;
line([mean(CVaRinN2{2,1}) mean(CVaRinN2{2,1})], [0 max(frequencyCounts)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')
hold on;
line([mean(CVaRinN2{2,2}) mean(CVaRinN2{2,2})], [0 max(frequencyCounts)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')
hold off;
text(1.2*mean(VaRinN2{2,1}), max(frequencyCounts)/1.5, 'VaR99%')
text(0.7*mean(CVaRinN2{2,1}), max(frequencyCounts)/1.3, 'CVaR99%')
text(1.1*mean(VaRinN2{2,2}), max(frequencyCounts)/1.4, 'VaR99.9%')
text(0.7*mean(CVaRinN2{2,2}), max(frequencyCounts)/1.1, 'CVaR99.9%')
title('MC2 loss distribution, portfolio 2(normal) - VaR and CVaR values')
xlabel('Loss $')
ylabel('Frequency')
