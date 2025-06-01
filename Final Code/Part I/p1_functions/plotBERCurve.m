function hFig = plotBERCurve(ber_data, params)
    EbN0_domain_dB = params.simulation.EbN0_domain_dB;
    ModType = params.modulation.ModulationType;
    ModOrder = params.modulation.ModulationOrder;

    if strcmpi(ModType, 'qam')
        BER_theoretical = berawgn(EbN0_domain_dB, 'qam', ModOrder, 'nondiff');
        theoryDisplayName = sprintf('Theoretical $%d$-QAM', ModOrder);
    elseif strcmpi(ModType, 'pam')
        BER_theoretical = berawgn(EbN0_domain_dB, 'pam', ModOrder, 'nondiff');
        if ModOrder == 2
            theoryDisplayName = 'Theoretical BPSK ($2$-PAM)';
        else
            theoryDisplayName = sprintf('Theoretical $%d$-PAM', ModOrder);
        end
    end

    hFig = figure('Name', 'BER Curve', 'NumberTitle', 'off');
    ax = gca; 

    semilogy(ax, EbN0_domain_dB, BER_theoretical, '-o', 'DisplayName', theoryDisplayName, 'LineWidth', 1);
    hold(ax, 'on');
    semilogy(ax, EbN0_domain_dB, ber_data, '-^', 'DisplayName', 'Simulated', 'LineWidth', 1);
    
    grid(ax, 'on');
    
    xlabel(ax, '$E_b/N_0$ (dB)', 'Interpreter', 'latex', 'FontSize', 30);
    ylabel(ax, 'Bit Error Rate (BER)', 'Interpreter', 'latex', 'FontSize', 30);
    
    % titleStr = sprintf('BER Performance for $%d$-%s (OSF=$%d$)', ModOrder, upper(ModType), params.sampling.OversamplingFactor);
    titleStr = sprintf('BER Performance for $%d$-%s', ModOrder, upper(ModType));
    title(ax, titleStr, 'Interpreter', 'latex', 'FontSize', 30);
    
    legend(ax, 'show', 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 20);
    
    ylim(ax, [1e-5 1]);
    xlim(ax, [min(EbN0_domain_dB) max(EbN0_domain_dB)]);
    
    ax.FontSize = 30;
    ax.TickLabelInterpreter = 'latex'; 
    
    hold(ax, 'off');
end