figure('name', 'Cumulative Acoustic energy vs Stress','Position',...
    [60,60,1400,700])
hold on
yyaxis left
[~, index] = max(AllValues.PARA1);
xline(AllValues.HitTime(index));
plot(time2/Resolution, SpreadEner);
% FOR POLYNOMIAL % 
% p = polyfit(time2/Resolution,SpreadEner,20);
% pol = polyval(p, time2/Resolution);
% plot(time2/Resolution, pol,'--');
title('Cumulative Acoustic energy vs Stress');
xlim([0 TimeEnd]);
%ylim([0 max(SpreadEner)+2]);
xlabel('Time [s]');
ylabel('Energy');
plot(time2/Resolution, SpreadEner1,'--');
plot(time2/Resolution, SpreadEner2,'--');
plot(time2/Resolution, SpreadEner3,'--');
xline(PullStop);
legend("Total", "0-200kHz", "200-400kHz", ">400kHz",...
    'location','south outside');
yyaxis right
% plot(AllValues.HitTime, AllValues.PARA1);
plot(CSVDataOffs.Fun_Time, CSVDataOffs.Fun_TensileStress);
hold off