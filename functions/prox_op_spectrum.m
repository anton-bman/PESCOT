function [alpha,dual0,dual1] = prox_op_spectrum(C,Phi,gradient,epsilon, zeta, eta, gamma,max_iter,dual0,dual1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   zeta: OT regularization parameter
%   eta: OT sparsity regularization parameter
%   epsilon: entropic regularization parameter

    r = Phi;

    r_filt = r>0;
    if sum(r_filt) == 0
        alpha = r*0;
        return
    end

    r_filtered = r(r_filt);
    C_filtered = C(r_filt,:);
    gradient_filtered = gradient(r_filt);

    
    dual1_full = zeros(size(dual1));
    dual0_full = -inf*ones(size(dual0));

    dual1 = dual1(r_filt,:);
    dual0 = dual0(r_filt);


    rel_error_0 = zeros(max_iter,1);
    rel_error_1 = zeros(max_iter,1);
    
    %vectors of 1
    ones_1 = ones(size(dual1,2),1);

    dual1_new = zeros(size(dual1));
    for iter = 1:max_iter
        %%%% ---- dual0 %%%%
        
        logxi = stable_log_sum_exp(-C_filtered/epsilon + dual1/(epsilon*zeta*gamma),ones_1);
        dual0_new = (log(r_filtered)-logxi -gamma*gradient_filtered )*(zeta*epsilon*gamma)/(1 + zeta*epsilon*gamma);



        %%%% ---- dual1 %%%%
        replicated_dual0_new = repmat(dual0_new, 1, size(C_filtered, 2));
        unsorted_mat = 1/(zeta*gamma*epsilon)*replicated_dual0_new - C_filtered/epsilon;
        sorted_transformed_matrix = sort(1/(zeta*gamma*epsilon)*replicated_dual0_new - C_filtered/epsilon, 1, 'ascend');
        
        for idx = 1:size(dual1_new,2)
            dual1_new(:,idx) = -(zeta*gamma*epsilon)*min_logsumexp_l1_ball_faster(sorted_transformed_matrix(:, idx),eta/epsilon, unsorted_mat(:, idx));
        end


        % Stopping criterion
        rel_error_0(iter) = norm(dual0-dual0_new)/norm(dual0_new);
        rel_error_1(iter) = norm(dual1-dual1_new, 'fro')/norm(dual1_new,'fro');

        %saving the new dual variables
        dual0 = dual0_new;
        dual1 = dual1_new;

        
        if iter > 1
            if rel_error_0(iter)<1e-3 && rel_error_1(iter)<1e-3
                break;
            end
        end
        if sum(isnan(dual1_new)) > 0
            fprintf('Numerical instability \n')
            return
        end
        if sum(isnan(rel_error_0(iter))) > 0
            fprintf('Numerical instability \n')
            alpha = r*0;
            break
        end
    end
    
    rhat = zeros(size(r));
    rhat(r_filt) = r(r_filt).*exp(-gamma*gradient_filtered-dual0);
    alpha = rhat;
    if rhat < 0
        print("r-hat less than 0, error")
        return;
    end

    dual1_full(r_filt,:) = dual1;
    dual0_full(r_filt) = dual0;

    dual1 = dual1_full;
    dual0 = dual0_full;
end