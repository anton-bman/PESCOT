function [alpha,dual0,dual1] = prox_operator_abs_l1(C,u,epsilon, zeta, eta, gamma,max_iter,dual0,dual1, beta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    omega = angle(u);
    r = abs(u);

    r_beta = r-beta*gamma;

    theta = zeta*gamma;

    rel_error_0 = zeros(max_iter,1);
    rel_error_1 = zeros(max_iter,1);
    
    ones_0 = ones(size(dual0,1),1);
    ones_1 = ones(size(dual1,2),1);

    dual0_new = zeros(size(dual0));
    dual1_new = zeros(size(dual1));
    for iter = 1:max_iter
        %%%% ---- dual0 %%%%


        xi = stable_log_sum_exp(-C/epsilon + dual1/(epsilon*theta),ones_1);
        dual_idx = log(r_beta) > xi;
        dual0_partial = max( r_beta(dual_idx) -epsilon*theta*wrightOmegaq(xi(dual_idx) - log(theta*epsilon) + r_beta(dual_idx)/(theta*epsilon)),0);
        dual0_new(~dual_idx) = 0;
        dual0_new(dual_idx) = dual0_partial;
        
        %%%% ---- dual1 %%%%
        replicated_dual0_new = repmat(dual0_new, 1, size(C, 2));
        unsorted_mat = 1/(theta*epsilon)*replicated_dual0_new - C/epsilon;
        sorted_transformed_matrix = sort(1/(theta*epsilon)*replicated_dual0_new - C/epsilon, 1, 'ascend');
        for idx = 1:size(dual1_new,2)
            dual1_new(:,idx) = -(theta*epsilon)*min_logsumexp_l1_ball_faster(sorted_transformed_matrix(:, idx),eta/epsilon, unsorted_mat(:, idx));
        end

        % Stopping criterion
        rel_error_0(iter) = norm(dual0-dual0_new)/norm(dual0_new);
        rel_error_1(iter) = norm(dual1-dual1_new, 'fro')/norm(dual1_new,'fro');

        % Saving the new dual variables
        dual0 = dual0_new;
        dual1 = dual1_new;

        if iter > 1
            if rel_error_0(iter)<1.5e-3 && rel_error_1(iter)<1.5e-3
                break;
            end
        end
        if sum(isnan(dual1_new)) > 0
            fprintf('Numerical instability \n')
            return
        end

    end
    
    rhat = max(r_beta - dual0,0);
    alpha = rhat.*exp(1j*omega);
    if rhat < 0
        print("r-hat less than 0, error")
        return;
    end
end