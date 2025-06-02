function [symb_rx_corrected_down, time_errors] = gardner(signal_rx, kappa, OSF)
   % Returns the downsampled signal. No need to downsample again in the
   % main code

    N_sym = floor(length(signal_rx) / OSF); 
    symb_rx_corrected_down = zeros(1, N_sym);
    time_errors = zeros(1, N_sym);

    symb_rx_corrected_down(1) = signal_rx(1);
    time_errors(1) = 0;
    
    for n = 1:N_sym - 1
        idx_range = (n - 1) * OSF + 1 : (n * OSF + 2);
        x = idx_range;
        symbols = signal_rx(x);

        eps_now = time_errors(n);
        idx_now = n * OSF + 1 - eps_now;
        idx_half = n * OSF + 1 - OSF/2 - eps_now;

        y_now = interp1(x, symbols, idx_now, 'linear', 0);
        y_prev = symb_rx_corrected_down(n);
        y_half = interp1(x, symbols, idx_half, 'linear', 0);

        eps_new = eps_now + kappa * real(y_half * (conj(y_now) - conj(y_prev)));
        symb_rx_corrected_down(n+1) = y_now;
        time_errors(n+1) = eps_new;
    end

    symb_rx_corrected_down = symb_rx_corrected_down.';
end