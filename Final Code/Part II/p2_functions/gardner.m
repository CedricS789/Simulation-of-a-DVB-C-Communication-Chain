function [symb_adjusted, time_errors] = gardner(y)
    % Gardner timing recovery for signal y with 2 samples per symbol

    k = 0.001;  % Step size
    N_sym = floor(length(y) / 2);  % Number of symbols (2 samples per symbol)

    time_errors = zeros(1, N_sym-1);
    symb_adjusted = zeros(1, N_sym);

    eps_now = 0;
    eps_prev = 0;

    for n = 2:N_sym
        % Compute sample positions
        idx_prev = 2*(n-2) + 1 - eps_prev;
        idx_half = 2*(n-1) - eps_now;
        idx_now  = 2*(n-1) + 1 - eps_now;

        % Interpolated values
        y_prev = interpolate(y, idx_prev);
        y_half = interpolate(y, idx_half);
        y_now  = interpolate(y, idx_now);

        % Gardner error and update
        error = real(y_half * (conj(y_now) - conj(y_prev)));
        eps_new = eps_now - k * error;

        % Store outputs
        time_errors(n-1) = eps_new;
        symb_adjusted(n-1) = y_now;

        % Update timing
        eps_prev = eps_now;
        eps_now = eps_new;
    end
end

