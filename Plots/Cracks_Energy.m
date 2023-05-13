cracks = [126, 122, 44, 6, 5, 0, 0];
p_cracks = [87, 79, 13, 1, 3, 0, 0];
NT_cracks = 258; 
NT_p_cracks = 87;

energy_titles = ["Pullstop Energy", "Pullstop + 10sec", "Total Energy"];
energy_stop = ["Var3", "Var4", "Var2"];

for t = energy_titles
    i = find(t == energy_titles);
    figure('name', plus("Matrix Cracks vs ", t),...
        'Position',[60,60,1400,700])
    hold on
    plot(SV.EnergyTable.energy_stop(i), cracks);
    title(append(plus("Matrix Cracks vs ", t),...
        captext));
    xlabel('Energy [aJ]');
    ylabel('No. of cracks');
    legend(Legendtext,'location','east outside');
    grid on
    set(gca,'FontSize',14)
    hold off
end