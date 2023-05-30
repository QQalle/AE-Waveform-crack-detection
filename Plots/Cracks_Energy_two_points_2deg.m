% Instructions:
% Enter all crack data below
% But the data of the specimens without tabs last
% Enter how many specimens didn't have tabs

% cracks = [126, 122, 44, 6, 5, 0, 0];
% p_cracks = [87, 79, 13, 1, 3, 0, 0];
cracks = [106.5, 100.5, 28.5, 3.5, 4, 0, 0];

point_one = 4;
point_two = 2;

no_tabs = 0;

energy_data = table2array(SV.EnergyTable(1:(end-no_tabs), 2:4));
energy_data_notabs = table2array(SV.EnergyTable((end-no_tabs+1):end, 2:4));

cracks_data = cracks(:, 1:(end-no_tabs));
cracks_data_notabs = cracks(:, (end-no_tabs+1):end);

max_stress = cellfun(@max, SV.Fun_TensileStress);

% energy_titles = ["Pullstop Energy", "Pullstop + 10sec", "Total Energy"];
% energy_stop_col = [2, 3, 1];

energy_titles = ["Pullstop Energy"];
energy_stop_col = [2];


for t = energy_titles

    i = find(t == energy_titles);
    e = energy_data(:, energy_stop_col(i))';
    e_nt = energy_data_notabs(:, energy_stop_col(i))';

    figure('name', append("Matrix Cracks vs ", t, ", using two points of data"),...
        'Position',[60,60,1400,700]);
    hold on

    lf_mc_en = polyfit([e(point_one) e(point_two)],...
        [cracks_data(point_one) cracks_data(point_two)], 1);
    % lcp = polyfit(x, p_cracks, 1);
    lff_mc_en = polyval(lf_mc_en, [min(e) max(e)]);
    % fp = polyval(lcp, x);

    scatter(e, cracks_data, 50, "filled", "Color", "#0072BD");

    scatter([e(point_two) e(point_one)],...
        [cracks_data(point_two) cracks_data(point_one)],...
        100, [0.4660 0.6740 0.1880], "o", "LineWidth", 1.5);

    plot([min(e) max(e)], lff_mc_en, "LineStyle", "--", "Color", "#0072BD", ...
        "LineWidth", 2);

    title(append(plus("Matrix Cracks vs ", t),...
        ", using two points of data", captext));
    xlabel('Energy [aJ]');
    ylabel('Average amount of Matrix Cracks per side, using two points of data');
    legend(["Real matrix crack data", "Selected Points" "Line"], ...
        'location','southoutside');
    grid on
    set(gca,'FontSize',14)
    hold off

    %% Energy vs Stress
    figure('name', plus(t, " vs Max Stress"),...
        'Position',[60,60,1400,700]);
    hold on

    %ES_range = 200:1:max(max_stress);
    st_range = min(max_stress):1:max(max_stress);

    lf_en_st = polyfit(max_stress, e, 2);
    lff_en_st = polyval(lf_en_st, st_range);

    scatter(max_stress, e, 50, "filled", "Color", "#0072BD");
    plot(st_range, lff_en_st, "LineStyle", "--",...
        "Color", "#0072BD", "LineWidth", 2);

    title(append(plus(t, " vs Max Stress"),...
        captext));
    xlabel('Max Stress [MPa]');
    ylabel('Accumulated Energy [aJ]');
    legend(["Accumulated Energy", "Polynomial fit (2 deg.)"], ...
        'location','southoutside');
    grid on
    set(gca,'FontSize',14)
    hold off

    %% Matrix Cracks vs Stress
    figure('name', "Matrix Cracks vs Max Stress, using two points of data", 'Position',[60,60,1400,700]);
    hold on

    max_stress_range = min(max_stress):1:max(max_stress);

    lff_mc_st = lf_mc_en(:,1) * (lf_en_st(:,1) * max_stress_range.^2 + lf_en_st(:,2) * max_stress_range + lf_en_st(:,3)) + lf_mc_en(:,2);

    %lf_mc_st_real = polyfit(max_stress, cracks, 2);
    %lff_mc_st_real = polyval(lf_mc_st_real, st_range);

    scatter(max_stress, cracks, 50, "filled", "Color", "#0072BD");

    plot(max_stress_range, lff_mc_st, "LineStyle", "--",...
        "Color", "#77AC30", "LineWidth", 2);

    title(append("Matrix Cracks vs Max Stress, using two points of data", captext));
    xlabel('Max Stress [MPa]');
    ylabel('Average amount of Matrix Cracks per side');
    legend(["Real matrix crack data", "Expected number off cracks (calculated using matrix crack data from two specimen)"], ...
        'location','southoutside');
    grid on
    set(gca,'FontSize',14)
    hold off

end





