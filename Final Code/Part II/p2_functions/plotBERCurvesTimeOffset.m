function hFig = plotBERCurvesTimeOffset(all_ber_data, params, time_offset_norm_values)
    EbN0_domain_dB = params.simulation.EbN0_domain_dB;
    ModType = params.modulation.ModulationType;
    ModOrder = params.modulation.ModulationOrder;

    if strcmpi(ModType, 'qam')
        BER_theoretical = berawgn(EbN0_domain_dB, 'qam', ModOrder, 'nondiff');
        theory_name = sprintf('Theoretical');
    elseif strcmpi(ModType, 'pam')
        BER_theoretical = berawgn(EbN0_domain_dB, 'pam', ModOrder, 'nondiff');
        if ModOrder == 2
            theory_name = 'Theoretical';
        else
            theory_name = sprintf('Theoretical');
        end
    end

    hFig = figure('Name', 'BER Curves vs. Sample Time Offset', 'NumberTitle', 'off');
    ax = gca;

    semilogy(ax, EbN0_domain_dB, BER_theoretical, 'DisplayName', theory_name, 'LineWidth', 1);
    hold(ax, 'on');
    
    markers = {'o', 's', 'd', '^', 'v', 'p', 'h', 'x', '+'};
    num_markers = length(markers);

    for i = 1:length(time_offset_norm_values)
        current_ber_data = all_ber_data(:, i);
        current_offset_val = time_offset_norm_values(i);
        simulation_name = sprintf('$t_0$ = %.2f$T_{symb}$', current_offset_val); 
        
        marker_style = markers{mod(i-1, num_markers) + 1}; 
        h_line_sim = semilogy(ax, EbN0_domain_dB, current_ber_data, 'Marker', marker_style, 'DisplayName', simulation_name, 'LineWidth', 1, 'MarkerSize', 7);
        set(h_line_sim, 'MarkerFaceColor', get(h_line_sim, 'Color'));
    end
    
    grid(ax, 'on');
    xlabel(ax, '$E_b/N_0$ (dB)', 'Interpreter', 'latex', 'FontSize', 30);
    ylabel(ax, 'Bit Error Rate (BER)', 'Interpreter', 'latex', 'FontSize', 30);
    titleStr = sprintf('BER vs. $E_b/N_0$ ($%d$-%s)', ModOrder, upper(ModType));
    title(ax, titleStr, 'Interpreter', 'latex', 'FontSize', 30);
    legend(ax, 'show', 'Location', 'southwest', 'Interpreter', 'latex', 'FontSize', 20);
    
    ylim(ax, [1e-5 1]);
    xlim(ax, [min(EbN0_domain_dB) max(EbN0_domain_dB)]);
    
    set(ax, 'YScale', 'log');
    ax.FontSize = 30;
    ax.TickLabelInterpreter = 'latex';
    ax.GridAlphaMode = 'manual'; 
    ax.GridAlpha = 0.6; 
    hold(ax, 'off');
end