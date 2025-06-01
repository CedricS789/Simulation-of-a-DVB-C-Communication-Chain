function hFig = plotConstellation_Tx_Rx(ModulationOrder, ModulationType, symb_tx, symb_rx)
    modStr = sprintf('%d-%s', ModulationOrder, upper(ModulationType));
    
    hFig = figure('Name', [modStr ' Constellation'], 'NumberTitle', 'off');
    set(hFig, 'DefaultTextInterpreter', 'latex', 'DefaultLegendInterpreter', 'latex', 'DefaultAxesTickLabelInterpreter', 'latex');
    
    ax = gca;
    ax.FontSize = 30;
    hold(ax, 'on');
    
    txColor = [0, 0.4470, 0.7410];
    rxColor = [0.8500, 0.3250, 0.0980];
    
    scatter(ax, real(symb_rx), imag(symb_rx), 60, 'x', 'MarkerEdgeColor', rxColor, 'LineWidth', 1.5, 'DisplayName', 'Received Symbols');
    scatter(ax, real(symb_tx), imag(symb_tx), 300, 'o', 'filled', 'MarkerFaceColor', txColor, 'MarkerEdgeColor', txColor, 'DisplayName', 'Transmitted Symbols');

   
    axis(ax, 'equal');
    xlim(ax, [-1.5, 1.5]);
    ylim(ax, [-1.5, 1.5]);

    xlabel(ax, '$\mathrm{In-Phase}$', 'FontSize', 30);
    ylabel(ax, '$\mathrm{Quadrature}$', 'FontSize', 30);
    % title(ax, sprintf('$%s$: Transmitted vs. Received Symbols', modStr), 'FontSize', 30);
    
    grid(ax, 'on');
    grid(ax, 'minor');
    lgd = legend(ax, 'show', 'Location', 'best');
    lgd.FontSize = 30;
    box(ax, 'on');
    hold(ax, 'off');
end