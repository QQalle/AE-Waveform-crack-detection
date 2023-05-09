figure('name', 'Hitcounter','Position',[60,60,1400,700])
disp("Spreading hits...")
hold on
%time = (0:ceil(max(HitTimeList)));
time = (0:TimeEnd);
SpreadHits = zeros(1,length(time));
for i = 1 : length(HitTimeList)
    SpreadHits(ceil(HitTimeList(i)):end) = ...
        SpreadHits(ceil(HitTimeList(i)):end)+1;
end
plot(time, SpreadHits);
title('Hits over time');
xlim([0 max(HitTimeList)]);
ylim([0 max(SpreadHits)+2]);
xlabel('Time [s]');
ylabel('Hits');
xline(PullStop);
DBSpreadHits = zeros(1,length(time)); %Debondings
if istable(Debondings) %If true = there are matrixcracks
    for i = 1 : length(Debondings.HitIndex)
        DBSpreadHits(ceil(Debondings.HitTime(i)):end) = ...
            DBSpreadHits(ceil(Debondings.HitTime(i)):end)...
            + 1;
    end
end
if istable(Debondings) %If true = there are debondings
    plot(time, DBSpreadHits,'g');
end
MCSpreadHits = zeros(1,length(time)); %Matrix cracks
if istable(Matrixcracks) %If true = there are matrixcracks
    for i = 1 : length(Matrixcracks.HitIndex)
        MCSpreadHits(ceil(Matrixcracks.HitTime(i)):end) = ...
            MCSpreadHits(ceil(Matrixcracks.HitTime(i)):end)...
            + 1;
    end
end
if istable(Matrixcracks) %If true = there are matrixcracks
    plot(time, MCSpreadHits,'r');
end
legend("Hits","Debondings","Matrix cracks",'location','south outside');
hold off