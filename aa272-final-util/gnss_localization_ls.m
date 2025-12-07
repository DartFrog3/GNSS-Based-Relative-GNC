function [x_est, b_est, num_iter, pdop] = gnss_localization_ls( ...
    x0, b0, x_sv_m, y_sv_m, z_sv_m, corr_pr_m, ...
    max_iterations, tol, print_solutions)
    % Newton-Raphson LS Solver as in hw

    if nargin < 8 || isempty(max_iterations)
        max_iterations = 20;
    end
    if nargin < 9 || isempty(tol)
        tol = 1.0; % [m]
    end
    if nargin < 10 || isempty(print_solutions)
        print_solutions = false;
    end

    x_sv_m = x_sv_m(:);
    y_sv_m = y_sv_m(:);
    z_sv_m = z_sv_m(:);
    corr_pr_m = corr_pr_m(:);

    x_est = x0(:);
    b_est = b0;

    num_iter = 0;

    if print_solutions
        fprintf("Beginning solver...\n");
    end

    while num_iter < max_iterations
        if print_solutions
            fprintf("Iteration %d:\n", num_iter);
        end

        G = get_geometry_matrix(x_est, x_sv_m, y_sv_m, z_sv_m);

        rho_th = get_theoretical_pseudoranges(x_est, b_est, x_sv_m, y_sv_m, z_sv_m);

        residuals = corr_pr_m - rho_th;
        d_est = pinv(G) * residuals;

        % Update state
        x_est = x_est + d_est(1:3);
        b_est = b_est + d_est(4);

        num_iter = num_iter + 1;

        if norm(d_est(1:3)) < tol && abs(d_est(4)) < tol
            if print_solutions
                fprintf("Reached solution within desired tolerance of %.3f m in %d steps.\n", ...
                        tol, num_iter);
            end
            pdop = get_PDOP(G);
            return;
        end
    end
    pdop = get_PDOP(G);
    if print_solutions
        fprintf("Maximum number of %d iterations reached.\n", max_iterations);
    end
end

%% Helpers

function G = get_geometry_matrix(x_est, x_sv_m, y_sv_m, z_sv_m)
    % Build geometry matrix

    dx = x_sv_m - x_est(1);
    dy = y_sv_m - x_est(2);
    dz = z_sv_m - x_est(3);

    los_vecs = [dx, dy, dz];
    los_norm = sqrt(sum(los_vecs.^2, 2));

    G = -los_vecs ./ los_norm;

    % clock column (all ones)
    G = [G, ones(size(G,1),1)];
end

function rho_th = get_theoretical_pseudoranges(x_est, b_est, x_sv_m, y_sv_m, z_sv_m)
    % theoretical geometric pr helper

    dx = x_sv_m - x_est(1);
    dy = y_sv_m - x_est(2);
    dz = z_sv_m - x_est(3);

    ranges = sqrt(dx.^2 + dy.^2 + dz.^2);
    rho_th = ranges + b_est;
end

function PDOP = get_PDOP(G)
    H = inv(transpose(G) * G);
    H_P = H(1:3, 1:3);
    PDOP = sum( diag(H_P).^2 );
end
