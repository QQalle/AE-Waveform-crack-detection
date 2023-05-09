figure('name', 'Energy vs Time','Position',[60,60,1400,700])
hold on
scatter(HitTimeList, HAFImpEnerList, 60, '.');
title('Energy vs time');
% ylim([0 3*10^5]);
xlabel('Time [s]');
ylabel("Energy [aJ]");
refl = refline(0,MCminEner); %minimum duration
refl.Color = 'r';
if istable(Matrixcracks) %If true = there are matrixcracks
    plot(Matrixcracks.HitTime, Matrixcracks.Energy, 'o');
end
plot(CheckVariableTable.HitTime, CheckVariableTable.Energy,'*');
hold off