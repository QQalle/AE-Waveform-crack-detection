% Instructions:
% Enter all crack data below
% But the data of the specimens without tabs last
% Enter how many specimens didn't have tabs

cracks = [126, 122, 44, 6, 5, 0, 0, 258];
% p_cracks = [87, 79, 13, 1, 3, 0, 0, 87];

no_tabs = 3;

energy_data = table2array(SV.EnergyTable(1:(end-no_tabs), 2:4));
energy_data_notabs = table2array(SV.EnergyTable((end-no_tabs+1):end, 2:4));

cracks_data = cracks(:, 1:(end-no_tabs));
cracks_data_notabs = cracks(:, (end-no_tabs+1):end);

energy_titles = ["Pullstop Energy", "Pullstop + 10sec", "Total Energy"];
energy_stop_col = [2, 3, 1];

for t = energy_titles

    i = find(t == energy_titles);
    e = energy_data(:, energy_stop_col(i))';
    e_nt = energy_data_notabs(:, energy_stop_col(i))';

    figure('name', plus("Matrix Cracks vs ", t),...
        'Position',[60,60,1400,700]);
    hold on

    lc = polyfit(e, cracks_data, 1);
    % lcp = polyfit(x, p_cracks, 1);
    f = polyval(lc, [min(e) max(e)]);
    % fp = polyval(lcp, x);

    scatter(e, cracks_data, 50, "filled", "Color", "#0072BD");
    scatter(e_nt, cracks_data_notabs, 50, "filled", "Color", "#808080");
    plot([min(e) max(e)], f, "LineStyle", "--", "Color", "#0072BD");
    % scatter(x, p_cracks, 50, "filled", "b");
    % plot(x, fp, "Color", "b", "LineStyle","--");

    title(append(plus("Matrix Cracks vs ", t),...
        captext));
    xlabel('Energy [aJ]');
    ylabel('Amount of Matrix Cracks');
    legend(["With tabs", "Without tabs", "Linear fit (w/ tabs)"], ...
        'location','east outside');
    grid on
    set(gca,'FontSize',14)
    hold off
end