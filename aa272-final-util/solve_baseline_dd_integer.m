function delta_r_m = solve_baseline_dd_integer(r_chief_est, r_sats_m, gnss_c, gnss_d, k_idx)
    % DDCP relative LS solver

    if nargin < 5 || isempty(k_idx)
        k_idx = size(gnss_c.pseudorange_m, 1); % Default to last entry
    end

    lambda = gnss_c.wavelength_m; 
    pr_c_all = gnss_c.pseudorange_m(k_idx, :).';
    pr_d_all = gnss_d.pseudorange_m(k_idx, :).';
    
    cp_c_cycles = gnss_c.carrier_phase_cycles(k_idx, :).';
    cp_d_cycles = gnss_d.carrier_phase_cycles(k_idx, :).';
    cp_c_m = cp_c_cycles * lambda;
    cp_d_m = cp_d_cycles * lambda;

    valid_c = find(pr_c_all > 0);
    valid_d = find(pr_d_all > 0);
    common_sats = intersect(valid_c, valid_d);
    Ns = length(common_sats);
    
    % Feasibility Check
    if Ns < 4
        delta_r_m = [0;0;0]; 
        fixed_ambiguities = [];
        return; 
    end
    
    pos_sats = r_sats_m(:, common_sats);
    pr_c = pr_c_all(common_sats);
    pr_d = pr_d_all(common_sats);
    cp_c = cp_c_m(common_sats);
    cp_d = cp_d_m(common_sats);
    
    % take 1st as ref
    ref_idx = 1; 
    
    other_indices = setdiff(1:Ns, ref_idx);
    N_others = length(other_indices);
    num_states = 3 + N_others;
    x_est = zeros(num_states, 1);
    
    % DD_N init (DD_Phase - DD_Code) / lambda
    for k = 1:N_others
        i = other_indices(k);
        dd_pr = (pr_d(i) - pr_c(i)) - (pr_d(ref_idx) - pr_c(ref_idx));
        dd_cp = (cp_d(i) - cp_c(i)) - (cp_d(ref_idx) - cp_c(ref_idx));
        
        x_est(3+k) = (dd_cp - dd_pr) / lambda;
    end
    
    % Get float solution
    for iter = 1:8
        H = zeros(2*N_others, num_states);
        y = zeros(2*N_others, 1);
        
        pos_ref = pos_sats(:, ref_idx);
        vec_ref_c = pos_ref - r_chief_est;
        u_ref = vec_ref_c / norm(vec_ref_c);
        
        r_dep_guess = r_chief_est + x_est(1:3);
        
        dd_geom_ref = norm(pos_ref - r_dep_guess) - norm(pos_ref - r_chief_est);
        
        row_count = 1;
        for k = 1:N_others
            i = other_indices(k);
            pos_i = pos_sats(:, i);
            vec_i_c = pos_i - r_chief_est;
            u_i     = vec_i_c / norm(vec_i_c);
            
            dd_geom_i = (norm(pos_i - r_dep_guess) - norm(pos_i - r_chief_est)) - dd_geom_ref;
            
            dd_pr_obs = (pr_d(i) - pr_c(i)) - (pr_d(ref_idx) - pr_c(ref_idx));
            dd_cp_obs = (cp_d(i) - cp_c(i)) - (cp_d(ref_idx) - cp_c(ref_idx));
            
            geom_H = (u_ref' - u_i'); 
            
            y(row_count)      = dd_pr_obs - dd_geom_i;
            H(row_count, 1:3) = geom_H;

            n_float = x_est(3+k);
            y(row_count+N_others)      = dd_cp_obs - dd_geom_i - (n_float * lambda);
            H(row_count+N_others, 1:3) = geom_H;
            H(row_count+N_others, 3+k) = lambda;
            
            row_count = row_count + 1;
        end
        
        % Weight matrix
        w_pr = 1 / 2.0^2; % ~ 2m std
        w_cp = 1 / 0.05^2; % ~5cm std
        W = diag([ones(N_others,1)*w_pr; ones(N_others,1)*w_cp]);
        
        dx = (H'*W*H) \ (H'*W*y);
        x_est = x_est + dx;
        
        if norm(dx(1:3)) < 1e-4, break; end
    end
    
    % rounding IAR
    float_ambigs = x_est(4:end);
    fixed_ambigs = round(float_ambigs);
    
    % resolve.
    H_fix = zeros(2*N_others, 3);
    y_fix = zeros(2*N_others, 1);
    
    r_dep_guess = r_chief_est + x_est(1:3);
    dd_geom_ref = norm(pos_ref - r_dep_guess) - norm(pos_ref - r_chief_est);
    
    for k = 1:N_others
        i = other_indices(k);
        pos_i = pos_sats(:, i);
        
        vec_i_c = pos_i - r_chief_est;
        u_i     = vec_i_c / norm(vec_i_c);
        
        dd_geom_i = (norm(pos_i - r_dep_guess) - norm(pos_i - r_chief_est)) - dd_geom_ref;
        geom_H = (u_ref' - u_i');
        
        dd_pr_obs = (pr_d(i) - pr_c(i)) - (pr_d(ref_idx) - pr_c(ref_idx));
        dd_cp_obs = (cp_d(i) - cp_c(i)) - (cp_d(ref_idx) - cp_c(ref_idx));

        y_fix(k) = dd_pr_obs - dd_geom_i;
        H_fix(k,:) = geom_H;
        y_fix(k+N_others) = dd_cp_obs - dd_geom_i - (fixed_ambigs(k) * lambda);
        H_fix(k+N_others,:) = geom_H;
    end
    
    dx_fix = (H_fix'*W*H_fix) \ (H_fix'*W*y_fix);
    delta_r_m = x_est(1:3) + dx_fix;
end