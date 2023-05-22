figure('name','Rise time vs time','Position',[60,60,1400,700])
hold on
scatter(HitTimeList, HAFImpRiseList, '.');
title('Rise time vs time');
xlabel('Time [s]');
ylabel("Rise time [μs]");
if istable(Matrixcracks) %If true = there are matrixcracks
    plot(Matrixcracks.HitTime, Matrixcracks.RiseTime, 'o');
end
plot(CheckVariableTable.HitTime, CheckVariableTable.RiseTime,'*');
hold off