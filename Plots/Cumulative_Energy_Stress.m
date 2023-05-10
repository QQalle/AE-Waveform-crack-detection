figure('name', 'Cumulative Acoustic energy vs Stress','Position',...
    [60,60,1400,700])
hold on

yyaxis left
[~, index] = max(AllValues.PARA1);
p = plot(time2/Resolution, SpreadEner, ...
    time2/Resolution, SpreadEner1,...
    time2/Resolution, SpreadEner2,...
    time2/Resolution, SpreadEner3);
% xline(AllValues.HitTime(index));
title('Cumulative Acoustic energy vs Stress');
xlim([0 TimeEnd]);
xlabel('Time [s]');
ylabel('Energy [aJ]');
xline(PullStop);
set(gca,'YColor','k')
% ytl = get(gca, 'YTick'); 

yyaxis right
ylabel('Stress [MPa]');
% plot(AllValues.HitTime, AllValues.PARA1);
plot(CSVDataOffs.Fun_Time, CSVDataOffs.Fun_TensileStress);
% ytr = get(gca, 'YTick'); 
% ytrv = linspace(min(ytr), max(ytr), numel(ytl)); % Create New Right Tick Values Matching Number Of Left Ticks
% ytrc = compose('%.2f',ytrv); % Tick Label Cell Array
% set(gca, 'YTick',ytrv, 'YTickLabel',ytrc)   
legend("Total energy", "0-200kHz", "200-400kHz", ">400kHz",...
    "Pull stop","Stress",'location','south outside');
set(gca,'YColor','k')
grid on
hold off