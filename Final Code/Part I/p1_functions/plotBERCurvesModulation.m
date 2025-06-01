function hFig = plotBERCurvesModulation(all_ber_data, params, nbps_values)
    EbN0_domain_dB = params.simulation.EbN0_domain_dB;

    hFig = figure('Name', 'BER Curves vs. Modulation Order', 'NumberTitle', 'off');
    ax = gca;
    hold(ax, 'on');
    
    markers = {'o', 's', 'd', '^', 'v', 'p', 'h', 'x', '+'};
    num_markers = length(markers);

    for i = 1:length(nbps_values)
        current_nbps = nbps_values(i);
        ModOrder = 2^current_nbps;
        current_ber_data = all_ber_data(:, i);
        
        if current_nbps == 1 || mod(current_nbps, 2) ~= 0 
            current_actual_mod_type = 'pam';
        else
            current_actual_mod_type = 'qam';
        end
        
        if strcmpi(current_actual_mod_type, 'pam') && ModOrder == 2
            displayName = sprintf('BPSK');
        else
            displayName = sprintf('$%d$-%s', ModOrder, upper(current_actual_mod_type));
        end
        
        marker_style = markers{mod(i-1, num_markers) + 1}; 
        h_line = semilogy(ax, EbN0_domain_dB, current_ber_data, 'Marker', marker_style, 'DisplayName', displayName, 'LineWidth', 1, 'MarkerSize', 7);
        set(h_line, 'MarkerFaceColor', get(h_line, 'Color'));
    end
    
    grid(ax, 'on');
    xlabel(ax, '$E_b/N_0$ (dB)', 'Interpreter', 'latex', 'FontSize', 30);
    ylabel(ax, 'Bit Error Rate', 'Interpreter', 'latex', 'FontSize', 30);
    title(ax, 'BER vs. $E_b/N_0$', 'Interpreter', 'latex', 'FontSize', 30); 
    legend(ax, 'show', 'Location', 'southwest', 'Interpreter', 'latex', 'FontSize', 24);
    
    ylim(ax, [1e-5 1]);
    xlim(ax, [min(EbN0_domain_dB) max(EbN0_domain_dB)]);
    
    set(ax, 'YScale', 'log');
    ax.FontSize = 30;
    ax.TickLabelInterpreter = 'latex';
    hold(ax, 'off');
end