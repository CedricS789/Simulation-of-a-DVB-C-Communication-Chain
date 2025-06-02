function hFig = plotBERCurvesCFO(all_ber_data, params, ppm_values)
    EbN0_domain_dB = params.simulation.EbN0_domain_dB;
    ModType = params.modulation.ModulationType;
    ModOrder = params.modulation.ModulationOrder;

    if strcmpi(ModType, 'qam')
        BER_theoretical = berawgn(EbN0_domain_dB, 'qam', ModOrder, 'nondiff');
        theory_name = sprintf('Theoretical');
    elseif strcmpi(ModType, 'pam')
        BER_theoretical = berawgn(EbN0_domain_dB, 'pam', ModOrder, 'nondiff');
        if ModOrder == 2
            theory_name = 'Theoretical BPSK ($2$-PAM) (No CFO)';
        else
            theory_name = sprintf('Theoretical $%d$-PAM (No CFO)', ModOrder);
        end
    end

    hFig = figure('Name', 'BER Curves vs CFO', 'NumberTitle', 'off');
    ax = gca;

    h_line = semilogy(ax, EbN0_domain_dB, BER_theoretical, 'DisplayName', theory_name, 'LineWidth', 1);
    set(h_line, 'MarkerFaceColor', get(h_line, 'Color'));

    hold(ax, 'on');
    
    markers = {'o', 's', 'd', '^', 'v', 'p', 'h', 'x', '+'};
    % markers = {'o'};
    num_markers = length(markers);
    num_cfo_values = length(ppm_values);

    for i = 1:num_cfo_values
        current_ber_data = all_ber_data(:, i);
        current_ppm = ppm_values(i);
        simulation_name = sprintf('CFO = $%g$ ppm', current_ppm);
        marker_style = [markers{mod(i-1, num_markers) + 1}]; 
        h_line = semilogy(ax, EbN0_domain_dB, current_ber_data, 'Marker', marker_style, 'DisplayName', simulation_name, 'LineWidth', 1, 'MarkerSize', 7);
        set(h_line, 'MarkerFaceColor', get(h_line, 'Color'));
    end
    
    grid(ax, 'on');
    xlabel(ax, '$E_b/N_0$ (dB)', 'Interpreter', 'latex', 'FontSize', 30);
    ylabel(ax, 'Bit Error Rate', 'Interpreter', 'latex', 'FontSize', 30);
    titleStr = sprintf('BER vs. $E_b/N_0$ ($%d$-%s)', ModOrder, upper(ModType));
    title(ax, titleStr, 'Interpreter', 'latex', 'FontSize', 30);
    legend(ax, 'show', 'Location', 'southwest', 'Interpreter', 'latex', 'FontSize', 30);
    
    ylim(ax, [1e-5 1]);
    xlim(ax, [min(EbN0_domain_dB) max(EbN0_domain_dB)]);
    
    ax.FontSize = 30;
    ax.TickLabelInterpreter = 'latex';
    ax.GridColor = [0.3 0.3 0.3];
    ax.GridAlpha = 0.7;
    hold(ax, 'off');
    box on;
end