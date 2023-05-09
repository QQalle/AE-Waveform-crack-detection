figure('name', 'Total Cumulative Acoustic Energy','Position',[60,60,1400,700])
hold on
plot(time2/Resolution, SpreadEner);
% FOR POLYNOMIAL % 
% p = polyfit(time2/Resolution,SpreadEner,20);
% pol = polyval(p, time2/Resolution);
% plot(time2/Resolution, pol,'--');
title('Cumulative Acoustic Energy');
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
hold off