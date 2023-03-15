function [ghat, y_residue_matrix, A_all_omega] = LeastSquares_2D(y_matrix, omega_est)

    [Nx, My, T] = size(y_matrix);
    NM_num = Nx * My;
    [K_est, ~] = size(omega_est);
    A_all_omega = zeros(NM_num, K_est);
    for k_idx = 1 : K_est
        omega_est_idx = omega_est(k_idx, :);
        xhat_vec_idx = exp((1j * ((0 : Nx - 1) - (Nx - 1) / 2)' * omega_est_idx(1))) / sqrt(Nx);
        yhat_vec_idx = exp((1j * ((0 : My - 1) - (My - 1) / 2)' * omega_est_idx(2))) / sqrt(My);
        atom_vec_est_idx = kron(yhat_vec_idx, xhat_vec_idx);
        A_all_omega(:, k_idx) = atom_vec_est_idx;
    end
    y_residue_matrix = zeros(Nx, My,T);
    ghat = zeros(K_est,T);
    for t = 1:T
        y_matrix_t = y_matrix(:,:,t);
        y_vector_t = y_matrix_t(:);
        ghat(:,t) = A_all_omega \ y_vector_t;
        y_residue_vec = y_vector_t - A_all_omega * ghat(:,t);
        y_residue_matrix(:,:,t) = reshape(y_residue_vec, Nx, My);
    end
    
   

end