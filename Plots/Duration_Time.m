figure('name', 'Duration vs Time','Position',[60,60,1400,700])
hold on
scatter(HitTimeList, HAFImpDurList, 60, '.');
title('Duration');
%xlim([0 TimeEnd]);
xlabel('Time [s]');
ylabel("Duration [μs]");
refl = refline(0,MCminDur); %minimum duration
refl.Color = 'r';
refl2 = refline(0,N/Fs*10^6);
refl2.Color = 'k';
if istable(Matrixcracks) %If true = there are matrixcracks
    plot(Matrixcracks.HitTime, Matrixcracks.Duration, 'o');
end
plot(CheckVariableTable.HitTime, CheckVariableTable.Duration,'*');
hold off