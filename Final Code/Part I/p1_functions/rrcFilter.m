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
           g_rrc = real(g_rrc).';
          
           % Normalize so that the combined pulse h = conv(g, g) has peak amplitude 1.
           % This ensures the signal amplitude at the sampling instant is correctly scaled.
           h_rc = conv(g_rrc, g_rrc);
           center_idx_rc = floor(length(h_rc)/2) + 1;   % Index of the peak for zero delay
           peak_val_rc = h_rc(center_idx_rc);
           g_rrc = g_rrc / sqrt(peak_val_rc);           % Corresponding g to the Normalized h      
           h_rc = conv(g_rrc, g_rrc);                   % Normalized h
    
    
            % Plots
            figure;
            ax1 = gca;
            plot(ax1, freq / 1e6, abs(H_RC), 'Color', [0 0.4470 0.7410], 'LineWidth', 2, 'Marker', 'o', 'MarkerFaceColor', [0 0.4470 0.7410], 'DisplayName', '$|H(f)|$');
            hold(ax1, 'on');
            plot(ax1, [f1_hz, f1_hz]/1e6, ylim(ax1), 'r--', 'LineWidth', 1.5, 'DisplayName', '$f_1$');
            plot(ax1, [-f1_hz, -f1_hz]/1e6, ylim(ax1), 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
            plot(ax1, [f2_hz, f2_hz]/1e6, ylim(ax1), 'g--', 'LineWidth', 1.5, 'DisplayName', '$f_2$');
            plot(ax1, [-f2_hz, -f2_hz]/1e6, ylim(ax1), 'g--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
            hold(ax1, 'off');
            grid(ax1, 'on');

            xlabel(ax1, 'Frequency (MHz)', 'Interpreter', 'latex', 'FontSize', 30);
            ylabel(ax1, 'Magnitude', 'Interpreter', 'latex', 'FontSize', 30);
            title(ax1, ['Raised Cosine Filter Discrete Frequency Points ($\beta$ = ', num2str(Beta), ', Taps = ', num2str(NumTaps), ')'], 'Interpreter', 'latex', 'FontSize', 30);
            legend(ax1, 'show', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 30);

            xlim(ax1, [-7, 7]);
            ylim(ax1, [-0.1, 2.2]*1e-7);
            
            ax1.GridLineWidth = 2;
            ax1.GridColor = [0.3 0.3 0.3];
            ax1.GridAlpha = 0.5;
            ax1.TickLabelInterpreter = 'latex';
            ax1.FontSize = 30;
            box(ax1, 'on');

            
            figure;
            ax2 = gca;
            Ts = 1 / Fs;
            t_rrc = (-(NumTaps-1)/2:(NumTaps-1)/2)*Ts;
            t_rc = (-(NumTaps-1):(NumTaps-1))*Ts;
            k_min = ceil((1 - center_idx_rc) / OSF);
            k_max = floor((length(h_rc) - center_idx_rc) / OSF);
            k_values = k_min:k_max;
            sampling_indices = center_idx_rc + k_values * OSF;
            h_Tsymb = h_rc(sampling_indices);
            k_Tsymb = (k_values * Tsymb);

            stem(ax2, k_Tsymb * 1e6, h_Tsymb, 'Color', [0 0.4470 0.7410], 'LineWidth', 4, 'MarkerFaceColor', [0 0.4470 0.7410], 'Marker', 'o', 'DisplayName', '$h(kT_{symb})$');
            hold(ax2, 'on');
            plot(ax2, t_rc * 1e6, h_rc, 'DisplayName', '$h(t)$', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);
            plot(ax2, t_rrc * 1e6, g_rrc, 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2, 'LineStyle','--', 'DisplayName', '$g(t)$');
            hold(ax2, 'off');
            grid(ax2, 'on');
            

            xlabel(ax2, 'Time ($\mu s$)', 'Interpreter', 'latex', 'FontSize', 30);
            ylabel(ax2, 'Amplitude', 'Interpreter', 'latex', 'FontSize', 30);
            title(ax2, 'Normalized RC and RRC Filter Impulse Responses', 'Interpreter', 'latex', 'FontSize', 30);
            legend(ax2, 'show', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 30);

            xlim(ax2, [-2, 2]);
            ylim(ax2, [-0.3 1.1]);
            
            ax2.GridLineWidth = 2;
            ax2.GridColor = [0.3 0.3 0.3];
            ax2.GridAlpha = 0.5;
            ax2.TickLabelInterpreter = 'latex';
            ax2.FontSize = 30;
            box(ax2, 'on');
        
end