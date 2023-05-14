cracks = [126, 122, 44, 6, 5, 0, 0];
p_cracks = [87, 79, 13, 1, 3, 0, 0];
NT_cracks = 258; 
NT_p_cracks = 87;

energy_data = table2array(SV.EnergyTable(1:(end-1), 2:4));
enrgy_data_notabs = table2array(SV.EnergyTable(end, 2:4));

energy_titles = ["Pullstop Energy", "Pullstop + 10sec", "Total Energy"];
energy_stop_col = [2, 3, 1];

for t = energy_titles

    i = find(t == energy_titles);
    x = energy_data(:, energy_stop_col(i))';

    figure('name', plus("Matrix Cracks vs ", t),...
        'Position',[60,60,1400,700])
    hold on

    lc = polyfit(x, cracks, 1);
    lcp = polyfit(x, p_cracks, 1);
    f = polyval(lc, [min(x) max(x)]);
    fp = polyval(lcp, x);

    scatter(x, cracks, 50, "filled", "Color", "#0072BD");
    scatter(enrgy_data_notabs(i), NT_cracks, 50, "filled", "Color", "#808080");
    plot([min(x) max(x)], f, "LineStyle", "--", "Color", "#0072BD");
%    scatter(x, p_cracks, 50, "filled", "b");
%    plot(x, fp, "Color", "b", "LineStyle","--");

    title(append(plus("Matrix Cracks vs ", t),...
        captext));
    xlabel('Energy [aJ]');
    ylabel('No. of cracks');
    legend(["Cracks w/ tabs", "Cracks w/o tabs", "Linear fit (w/ tabs)"], ...
        'location','east outside');
    grid on
    set(gca,'FontSize',14)
    hold off
end