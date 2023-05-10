figure('name', 'A-E vs Stress 400-infkHz','Position',...
    [60,60,1400,700])
hold on
for k = 1 : exp
%     yyaxis left
%     plot(SV.time2{k}/SV.Resolution(k), SV.SpreadEner{k});
    plot(SV.Fun_TensileStress{k}, SV.SpreadEner3{k});
%     yyaxis right
%     plot(SV.Fun_Time{k}, SV.Fun_TensileStress{k})
end
title('Cumulative Acoustic Energy vs Stress 400-infkHz');
% xlim([0 max(SV.Fun_Time)]);
% xlabel('Time [s]');
xlabel('Stress [MPa]');
ylabel('Energy [aJ]');
legend(num2str(round(SV.Table.MPa)));
hold off