function g_rrc = rrcFilter(Beta, SymRate, OSF, NumTaps)
           Fs = OSF * SymRate;
           Tsymb = 1 / SymRate;
           
           freq = linspace(-Fs/2, Fs/2, NumTaps);
           f1_hz = (1 - Beta) / (2 * Tsymb);
           f2_hz = (1 + Beta) / (2 * Tsymb);
           H_RC = zeros(1, NumTaps);
       
           % Region 1: Passband (|f| <= f1)
           pass_band_idx = (abs(freq) <= f1_hz);
           H_RC(pass_band_idx) = Tsymb;
       
           % Region 2: Roll-off band (f1 < |f| <= f2)
           roll_off_idx = (abs(freq) > f1_hz) & (abs(freq) <= f2_hz);
           H_RC(roll_off_idx) = (Tsymb / 2) * (1 + cos((pi * Tsymb / Beta) * (abs(freq((roll_off_idx))) - f1_hz)));
       
           % Calculate RRC frequency response: G(f) = sqrt(H(f))
           G_RRC = sqrt(H_RC);
       
           % IFFT to get time-domain impulse response, centered using fftshift/ifftshift
           g_rrc = fftshift(ifft(ifftshift(G_RRC)));
       
           % Ensure real output (due to potential numerical inaccuracies)
           g_rrc = real(g_rrc)';
          
           % Normalize so that the combined pulse h = conv(g, g) has peak amplitude 1.
           % This ensures the signal amplitude at the sampling instant is correctly scaled.
           h_rc = conv(g_rrc, g_rrc);
           center_idx_rc = floor(length(h_rc)/2) + 1;   % Index of the peak for zero delay
           peak_val_rc = h_rc(center_idx_rc);
           g_rrc = g_rrc / sqrt(peak_val_rc);           % Corresponding g to the Normalized h      
           h_rc = conv(g_rrc, g_rrc);                   % Normalized h
    
    
           % Plots
            % figure;
            % set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 30);
            % plot(freq / 1e6, abs(H), 'Color', [0 0.4470 0.7410], 'LineWidth', 2, 'Marker', 'o', 'MarkerFaceColor', [0 0.4470 0.7410], 'DisplayName', '$|H(f)|$');
            % hold on;
            % plot([f1_hz, f1_hz]/1e6, ylim, 'r--', 'LineWidth', 1.5, 'DisplayName', '$f_1$');
            % plot([-f1_hz, -f1_hz]/1e6, ylim, 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
            % plot([f2_hz, f2_hz]/1e6, ylim, 'g--', 'LineWidth', 1.5, 'DisplayName', '$f_2$');
            % plot([-f2_hz, -f2_hz]/1e6, ylim, 'g--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
            % hold off;
            % grid on;
            % set(gca, 'GridLineWidth', 2, 'GridColor', [0.3 0.3 0.3], 'GridAlpha', 0.5);
            % xlabel('Frequency (MHz)', 'Interpreter', 'latex', 'FontSize', 30);
            % ylabel('Magnitude', 'Interpreter', 'latex', 'FontSize', 30);
            % title(['Raised Cosine Filter Discrete Frequency Points ($\beta$ = ', num2str(Beta), ', Taps = ', num2str(NumTaps), ')'], 'Interpreter', 'latex', 'FontSize', 30);
            % legend('show', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 30);
            % % xlim(Fs/2*[-1 1] / 1e6);
            % xlim([-7, 7]);
            % ylim([-0.1, 2.2]*1e-7)
            % box on;
            % 
            % 
            % 
            % figure;
            % set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 30);
            % Ts = 1 / Fs;
            % t_rrc = (-(NumTaps-1)/2:(NumTaps-1)/2)*Ts;
            % t_rc = (-(NumTaps-1):(NumTaps-1))*Ts;
            % k_min = ceil((1 - center_idx_rc) / OSF);
            % k_max = floor((length(h) - center_idx_rc) / OSF);
            % k_values = k_min:k_max;
            % sampling_indices = center_idx_rc + k_values * OSF;
            % h_Tsymb = h(sampling_indices);
            % k_Tsymb = (k_values * Tsymb);
            % stem(k_Tsymb * 1e6, h_Tsymb, 'Color', [0 0.4470 0.7410], 'LineWidth', 4, 'MarkerFaceColor', [0 0.4470 0.7410], 'Marker', 'o', 'DisplayName', '$h(kT_{symb})$');
            % hold on;
            % plot(t_rc * 1e6, h, 'DisplayName', '$h(t)$', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);
            % hold on;
            % plot(t_rrc * 1e6, g, 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2, 'LineStyle','--', 'DisplayName', '$g(t)$');
            % hold off;
            % grid on;
            % set(gca, 'GridLineWidth', 2, 'GridColor', [0.3 0.3 0.3], 'GridAlpha', 0.5);
            % xlabel('Time ($\mu s$)', 'Interpreter', 'latex', 'FontSize', 30);
            % ylabel('Amplitude', 'Interpreter', 'latex', 'FontSize', 30);
            % title('Normalized RC and RRC Filter Impulse Responses', 'Interpreter', 'latex', 'FontSize', 30);
            % legend('show', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 30);
            % xlim([-5, 5]);
            % ylim([-0.3 1.1]);
            % box on;
        
end