close all
clear

experimentNo = [];
Experiments = ["2001","2002","2003","2004","2005",...
    "2006","2007","2008","2009"];
Experiments = ["2001","2002"];
SV = struct('experimentNo',experimentNo,'Experiments',Experiments);
for exp = 1 : length(SV.Experiments)
        %Variables to update
    experimentNo = SV.Experiments(exp);
    SV.experimentNo = experimentNo; %Variables to save
%     disp(experimentNo)
    if exp == 1
        StressThres = [];
        TEnerThres = [];
        SpreadEner = [];
        MPa = [];
        Table = table(MPa,StressThres,TEnerThres,SpreadEner);
        SV.Table = Table;
        SV.time2 = cell(1,length(SV.Experiments));
        SV.SpreadEner = cell(1,length(SV.Experiments));
        SV.SpreadEner1 = cell(1,length(SV.Experiments));
        SV.SpreadEner2 = cell(1,length(SV.Experiments));
        SV.SpreadEner3 = cell(1,length(SV.Experiments));
        SV.EnergyTable = table();
        SV.EnergyTable1 = table();
        SV.EnergyTable2 = table();
        SV.EnergyTable3 = table();
        SV.TimeEnd = [];
        SV.Fun_Time = cell(1,length(SV.Experiments));
        SV.Fun_TensileStress = cell(1,length(SV.Experiments));
        SV.Resolution = [];
        SV.HitTimeList = [];
        SV.SpreadEnerAPE = [];
        SV.SpreadEnerAPETime = [];
        SV.PullStop = [];
        SV.PullStopIndex = [];
        SV.PullStopIndexStress = [];
    end
    run Start.m
    SV.Table.MPa(exp) = max(CSVDataOffs.Fun_TensileStress);
    if isempty(StressThres) == 0
        SV.Table.StressThres(exp) = StressThres;
        SV.Table.TEnerThres(exp) = TEnerThres;
        SV.Table.SpreadEner(exp) = max(SpreadEner);
    end
    SV.time2{exp} = time2;
    SV.SpreadEner{exp} = SpreadEner;
    SV.SpreadEner1{exp} = SpreadEner1;
    SV.SpreadEner2{exp} = SpreadEner2;
    SV.SpreadEner3{exp} = SpreadEner3;
    SV.TimeEnd(exp) = TimeEnd;
    SV.Fun_Time{exp} = CSVDataOffs.Fun_Time;
    SV.Fun_TensileStress{exp} = CSVDataOffs.Fun_TensileStress;
    SV.Resolution(exp) = Resolution;
    SV.HitTimeList{exp} = HitTimeList;
    SV.SpreadEnerAPE{exp} = SpreadEnerAPE;
    SV.SpreadEnerAPETime{exp} = SpreadEnerAPETime;
    SV.PullStop(exp) = PullStop;
    SV.PullStopIndex(exp) = find(SV.HitTimeList{exp} >= PullStop,1);
    SV.PullStopIndexStress(exp) = find(SV.Fun_Time{exp} >= SV.PullStop(exp), 1);
    SV.EnergyTable{exp,1} = experimentNo;
    SV.EnergyTable{exp,2} = max(SV.SpreadEner{exp});
    SV.EnergyTable{exp,3} = SV.SpreadEner{exp}(SV.PullStopIndex(exp));
    SV.EnergyTable1{exp,1} = experimentNo;
    SV.EnergyTable1{exp,2} = max(SV.SpreadEner1{exp});
    SV.EnergyTable1{exp,3} = SV.SpreadEner1{exp}(SV.PullStopIndex(exp));
    SV.EnergyTable2{exp,1} = experimentNo;
    SV.EnergyTable2{exp,2} = max(SV.SpreadEner2{exp});
    SV.EnergyTable2{exp,3} = SV.SpreadEner2{exp}(SV.PullStopIndex(exp));
    SV.EnergyTable3{exp,1} = experimentNo;
    SV.EnergyTable3{exp,2} = max(SV.SpreadEner3{exp});
    SV.EnergyTable3{exp,3} = SV.SpreadEner3{exp}(SV.PullStopIndex(exp));
end
%%
Markers = ["o-","*-","x-","square-","diamond-","^-","v-","<-",...
    ">-","pentagram-","hexagram-"];
Markers = [Markers Markers Markers Markers];

figure('name', 'Cumulative Acoustic Energy vs Stress','Position',...
    [60,60,1400,700])
hold on
for k = 1 : exp
%     yyaxis left
%     plot(SV.time2{k}/SV.Resolution(k), SV.SpreadEner{k});
%     PullStopInd = find(SV.Fun_Time{k} >= SV.PullStop(k), 1);
    plot(SV.Fun_TensileStress{k}, SV.SpreadEner{k},Markers(k),...
        'MarkerIndices',SV.PullStopIndexStress(exp));
%     yyaxis right
%     plot(SV.Fun_Time{k}, SV.Fun_TensileStress{k})
end
title('Cumulative Acoustic Energy vs Stress');
% xlim([0 max(SV.Fun_Time)]);
% xlabel('Time [s]');
xlabel('Stress [MPa]');
ylabel('Energy [aJ]');
legend(num2str(round(SV.Table.MPa)));
hold off
%%
% Accumulated Energy after PullStop v. Time
figure('name', 'Proportion of total accumulated Energi after pullstop vs Time',...
    'Position',[60,60,1400,700])
hold on
for k = 1 : exp
    plot(SV.SpreadEnerAPETime{k}, 100*SV.SpreadEnerAPE{k});
end
title('Proportion of total accumulated energi After pullstop vs Time');
% xlim([0 max(SV.Fun_Time)]);
% xlabel('Time [s]');
xlabel('Time [sec]');
ylabel('Percentage of total Energy');
legend(plus(num2str(round(SV.Table.MPa)),"MPa"));
grid on
ytickformat('percentage')
hold off

% Stress v. Hits
figure('name','Stress vs Hits','Position',[60,60,1400,700])
hold on
for k = 1 : exp
    SpreadHits2 = zeros(1,length(SV.time2{k}));
    for i = 1 : length(SV.HitTimeList{k})
        SpreadHits2(ceil(SV.HitTimeList{k}(i))*SV.Resolution(k):end) = ...
            SpreadHits2(ceil(SV.HitTimeList{k}(i))*SV.Resolution(k):end)+1;
    end
    SmoothSpreadHits2 = smooth(SpreadHits2,100);
    plot(SV.Fun_TensileStress{k}, SmoothSpreadHits2);
end
title('Stress vs Hits');
xlabel('Stress [MPa]');
ylabel('Hits');
legend(plus(num2str(round(SV.Table.MPa)),"MPa"));
grid on
hold off

% Hits v. time (after pullstop)
figure('name','Hits vs Time','Position',[60,60,1400,700])
hold on
for k = 1 : exp
    SpreadHits2 = zeros(1,length(SV.time2{k}));
    for i = 1 : length(SV.HitTimeList{k})
        SpreadHits2(ceil(SV.HitTimeList{k}(i))*SV.Resolution(k):end) = ...
            SpreadHits2(ceil(SV.HitTimeList{k}(i))*SV.Resolution(k):end)+1;
    end
    SmoothSpreadHits2 = smooth(SpreadHits2,100);
    PullStopInd = find(SV.Fun_Time{k} >= SV.PullStop(k), 1);
    SmoothSpreadHits2origo = (SmoothSpreadHits2(PullStopInd:end)...
        - min(SmoothSpreadHits2(PullStopInd:end)))...
        ./ max(SmoothSpreadHits2);
    plot(SV.Fun_Time{k}(1:length(SmoothSpreadHits2origo))...
        , SmoothSpreadHits2origo);
end
title('Hits vs Time');
xlabel('Time [S]');
ylabel('Hits');
legend(plus(num2str(round(SV.Table.MPa)),"MPa"));
grid on
hold off

run Plots_cycle/A_E_Stress_0_200kHz.m
run Plots_cycle/A_E_Stress_200_400kHz.m
run Plots_cycle/A_E_Stress_400_infkHz.m

% SV.Table = sortrows(SV.Table);
