function hFig = plotConstellation_Tx_Rx(ModulationOrder, ModulationType, symb_tx, symb_rx)
    modStr = sprintf('%d-%s', ModulationOrder, upper(ModulationType));
    
    hFig = figure('Name', [modStr ' Constellation'], 'NumberTitle', 'off');
    set(hFig, 'DefaultTextInterpreter', 'latex', 'DefaultLegendInterpreter', 'latex', 'DefaultAxesTickLabelInterpreter', 'latex');
    
    ax = gca;
    ax.FontSize = 30;
    hold(ax, 'on');
    
    txColor = [0, 0.4470, 0.7410];
    rxColor = [0.8500, 0.3250, 0.0980];

    scatter(ax, real(symb_rx), imag(symb_rx), 60, 'x', 'MarkerEdgeColor', rxColor, 'LineWidth', 1.5, 'DisplayName', '$\mathrm{Received}$');
    scatter(ax, real(symb_tx), imag(symb_tx), 300, 'o', 'filled', 'MarkerFaceColor', txColor, 'MarkerEdgeColor', txColor, 'DisplayName', '$\mathrm{Transmitted}$');

    all_vals_for_max_calc = [real(symb_tx(:)); real(symb_rx(:)); imag(symb_tx(:)); imag(symb_rx(:))];
    max_coord = max([0; abs(all_vals_for_max_calc)], [], 'omitnan');
    
    lim_val_intermediate = max_coord * 1.15;
    lim_val = (lim_val_intermediate < 1e-5) * 1.0 + (lim_val_intermediate >= 1e-5) * lim_val_intermediate;

    plot(ax, [-lim_val, lim_val], [0 0], 'k--', 'LineWidth', 0.5, 'HandleVisibility','off');
    plot(ax, [0 0], [-lim_val, lim_val], 'k--', 'LineWidth', 0.5, 'HandleVisibility','off');

    axis(ax, 'equal');
    xlim(ax, [-lim_val, lim_val]);
    ylim(ax, [-lim_val, lim_val]);

    xlabel(ax, '$\mathrm{In-Phase} \; (\mathrm{I})$', 'FontSize', 30);
    ylabel(ax, '$\mathrm{Quadrature} \; (\mathrm{Q})$', 'FontSize', 30);
    title(ax, sprintf('$%s$: Transmitted vs. Received', modStr), 'FontSize', 30);
    
    grid(ax, 'on');
    grid(ax, 'minor');
    lgd = legend(ax, 'show', 'Location', 'best');
    lgd.FontSize = 20;
    box(ax, 'on');
    hold(ax, 'off');
end