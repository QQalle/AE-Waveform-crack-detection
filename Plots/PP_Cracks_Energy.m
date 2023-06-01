% Instructions:
% Enter all crack data below
% But the data of the specimens without tabs last
% Enter how many specimens didn't have tabs

% cracks = [126, 122, 44, 6, 5, 0, 0];
% p_cracks = [87, 79, 13, 1, 3, 0, 0];
cracks = [106.5, 100.5, 28.5, 3.5, 4, 0, 0];

no_tabs = 0;

energy_data = table2array(SV.EnergyTable(1:(end-no_tabs), 2:4));
energy_data_notabs = table2array(SV.EnergyTable((end-no_tabs+1):end, 2:4));

cracks_data = cracks(:, 1:(end-no_tabs));
cracks_data_notabs = cracks(:, (end-no_tabs+1):end);

max_stress = cellfun(@max, SV.Fun_TensileStress);

energy_titles = ["Pullstop Energy"];
energy_stop_col = [2];


for t = energy_titles

    i = find(t == energy_titles);
    e = energy_data(:, energy_stop_col(i))';
    e_nt = energy_data_notabs(:, energy_stop_col(i))';

    figure('name', "Matrix Cracks vs Acoustic Energy",...
        'Position',[60,60,1400/2,700/1.5]);
    hold on

    lf_mc_en = polyfit(e, cracks_data, 1);
    % lcp = polyfit(x, p_cracks, 1);
    lff_mc_en = polyval(lf_mc_en, [min(e) max(e)]);
    % fp = polyval(lcp, x);

    scatter(e, cracks_data, 50, "filled", "Color", "#0072BD");
    % scatter(e_nt, cracks_data_notabs, 50, "filled", "Color", "#808080");
    plot([min(e) max(e)], lff_mc_en, "LineStyle",...
        "--", "Color", "#0072BD", "LineWidth", 2);
    % scatter(x, p_cracks, 50, "filled", "b");
    % plot(x, fp, "Color", "b", "LineStyle","--");

    title("Matrix Cracks vs Acoustic Energy");
    xlabel('Energy [aJ]');
    ylabel('Number of Cracks');
    legend(["Data", "Linear fit"], ...
        'location','southoutside');
    grid on
    set(gca,'FontSize',14)
    hold off
end





