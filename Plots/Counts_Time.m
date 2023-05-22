figure('name', 'Counts vs Time','Position',[60,60,1400,700]);
hold on
scatter(HitTimeList, HAFImpCountList,'.');
title('Counts vs Time');
xlabel('Time [s]');
ylabel("Counts");
if istable(Matrixcracks) %If true = there are matrixcracks
    plot(Matrixcracks.HitTime, Matrixcracks.Counts, 'o');
end
plot(CheckVariableTable.HitTime, CheckVariableTable.Counts,'*');
hold off