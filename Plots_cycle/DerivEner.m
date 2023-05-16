figure('name', 'Energy-time derivative','Position',[60,60,1400,700])
disp("calculating Energy-time derivative...")
hold on

for k = 1 : exp
    PullStopInd = find(SV.Fun_TensileStress{k} >= SV.Table{k,1}, 1);
    plot(SV.Fun_TensileStress{k}, SV.DerivEner{k},Markers(k),...
        'MarkerIndices',PullStopInd,'MarkerSize',10);
%     disp(k);
%     disp(append('Peaks: ',num2str(EnerPeaksNo)));
%     disp(append('Peaks divided: ',num2str(EnerPeaksNo2)));
end
title('Energy Derivative vs Stress');
legend(Legendtext,'location','east outside');
xlabel('Stress [MPa]');
ylabel('Energy derivative');
set(gca,'FontSize',14)
hold off