function [y_residue_matrix, omega_hat, ghat] = RefineOne_2D(y_matrix, omega_est, ghat)
%UNTITLED2 此处显示有关此函数的摘要
%   date: 2022.1.11
    y_residue_matrix = y_matrix;
    obj_old = norm(y_residue_matrix(:));
    [Nx, My,T] = size(y_matrix);
    ant_idx_Nx = (0 : (Nx - 1))' - (Nx - 1) / 2;
    ant_idx_My = (0 : (My - 1))' - (My - 1) / 2;

    xhat_vec = exp((1j * ant_idx_Nx * omega_est(1))) / sqrt(Nx);
    yhat_vec = exp((1j * ant_idx_My * omega_est(2))) / sqrt(My);

%     atom_vec_est = kron(yhat_vec, xhat_vec);
%     atom_vec_est_matrix = reshape(atom_vec_est, Nx, My);
    y_matrix_add = zeros(Nx, My,T);
    x1diff_vec = 1j * ant_idx_Nx .* xhat_vec;
    y1diff_vec = 1j * ant_idx_My .* yhat_vec;

    % get the gradient of vector a
    ax1diff_vec = kron(yhat_vec, x1diff_vec);
    ay1diff_vec = kron(y1diff_vec, xhat_vec);

    da_over_domega = [ax1diff_vec, ay1diff_vec];
    Hessian_S = zeros(2, 2);
    grad_S = zeros(2, 1);
    for t = 1:T
        y_matrix_add(:,:,t) = y_matrix(:,:,t) + ghat(t) * xhat_vec*(yhat_vec.');
        y_matrix_add_temp = y_matrix_add(:,:,t);
        y_vector_add_temp = y_matrix_add_temp(:);
        grad_S = grad_S- 2 * real(conj(ghat(t)) * da_over_domega' * y_vector_add_temp);
        x2diff_vec = 1j * ant_idx_Nx .* x1diff_vec;
        y2diff_vec = 1j * ant_idx_My .* y1diff_vec;
        ax2diff_vec = kron(yhat_vec, x2diff_vec);
        ay2diff_vec = kron(y2diff_vec, xhat_vec);
        axy2diff_vec = kron(y1diff_vec, x1diff_vec);
        Hessian_S(1, 1) = Hessian_S(1, 1)- 2 * real(ghat(t)' * ax2diff_vec' * y_vector_add_temp);
        Hessian_S(2, 2) = Hessian_S(2, 2)- 2 * real(ghat(t)' * ay2diff_vec' * y_vector_add_temp);
        Hessian_S(1, 2) = Hessian_S(1, 2)- 2 * real(ghat(t)' * axy2diff_vec' * y_vector_add_temp);
    end
    Hessian_S(2, 1) = Hessian_S(1, 2);
%     y_matrix_add = y_matrix + ghat * atom_vec_est_matrix;
    

    Step = Hessian_S \ grad_S;
    if max(abs(Step)) < 2 * pi / min(Nx, My) / 4
        omega_hat_new = omega_est - Step';
        % sprintf('Newton is effective')
    else
        omega_hat_new = omega_est-sign(grad_S') * 2 * pi / max(Nx, My) / 100;
        % sprintf('Newton is ineffective and Grad is used')
    end

    xhat_vec_new = exp((1j * ant_idx_Nx * omega_hat_new(1))) / sqrt(Nx);
    yhat_vec_new = exp((1j * ant_idx_My * omega_hat_new(2))) / sqrt(My);
    ghat_new = zeros(1,T);
    for t = 1:T
        ghat_new(t) = xhat_vec_new'*y_matrix_add(:,:,t)*conj(yhat_vec_new);
        y_matrix(:,:,t) = y_matrix_add(:,:,t) - ghat_new(t) * xhat_vec_new*(yhat_vec_new.');
    end


%     atom_vec_est_new = kron(yhat_vec_new, xhat_vec_new);
%     ghat_new = atom_vec_est_new' * y_vector_add;
    
    obj_new = norm(y_matrix(:));
    if obj_new < obj_old
        y_residue_matrix = y_matrix;
        omega_hat = omega_hat_new;
        ghat = ghat_new;
        % sprintf('New is used')
    else
        y_residue_matrix = y_matrix;
        omega_hat = omega_est;
    end
end

