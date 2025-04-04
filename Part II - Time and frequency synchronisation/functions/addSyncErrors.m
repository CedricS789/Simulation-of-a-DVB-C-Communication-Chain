function signal_out = addSyncErrors(signal_in, delta_cfo_hz, phase_offset, Ts)
    % ---- signal_in: Input complex baseband signal vector (pulse-shaped). ---

    % ---- CFO Parameters ----
    delta_omega_offset = 2 * pi * delta_cfo_hz;               % Frequency offset in rad/s
    
    % --- Add CFO to the transmitted signal ---
    num_samples_tx  = length(signal_in);
    time_vector     = (0 : num_samples_tx - 1).' * Ts;          
    offset_signal   = exp(1j * delta_omega_offset * time_vector + phase_offset);
    signal_out      = signal_in .* offset_signal;              % Apply CFO to the transmitted signal

end