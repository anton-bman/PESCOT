function [alpha,dual0,dual1] = proximal_operator_square_l1(C,u,epsilon, zeta, eta, gamma,max_iter,dual0,dual1, beta)
% Proximal operator in the paper

    omega = angle(u);
    r = abs(u);

    r_beta = r-beta*gamma;
    r_beta_geq_0 = r-beta*gamma > 0;
    r_beta_positive = max(r_beta,0);

    theta = zeta*gamma;

    rel_error_0 = zeros(max_iter,1);
    rel_error_1 = zeros(max_iter,1);
    
    %vectors of 1
    ones_0 = ones(size(dual0,1),1);
    ones_1 = ones(size(dual1,2),1);

    dual1_new = zeros(size(dual1));
    for iter = 1:max_iter
        %%%% ---- dual0 %%%%
        xi = stable_log_sum_exp(-C/epsilon + dual1/(epsilon*theta),ones_1);
        % dual0_new = max(2*epsilon*theta*wrightOmegaq(1/2*log(r_beta.^2) -1/2*log(16*(epsilon*theta)^2) - ...
            % 1/2*stable_log_sum_exp(-C/epsilon + dual1/(epsilon*theta),ones_0) + 1/4/epsilon/theta) - 1/2,0);
            dual0_new = max(2*epsilon*theta*wrightOmegaq(1/2*log( ((r_beta).^2 - (  r_beta - r_beta_positive).^2) ) - 1/2*xi -log(4*theta*epsilon) +1/(4*theta*epsilon)  )-1/2,0);
        if dual0_new(~r_beta_geq_0) ~= 0
            error('dual variable wrong update')
        end
        
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
            if rel_error_0(iter)<0.001 && rel_error_1(iter)<0.001
                break;
            end
        end
        if sum(isnan(dual1_new)) > 0
            fprintf('Numerical instability \n')
            return
        end
    end
    
    rhat = r_beta_positive./(1+2*dual0);
    alpha = rhat.*exp(1j*omega);
    if rhat < 0
        print("r-hat less than 0, error")
        return;
    end
end