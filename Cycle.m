close all
clear

experimentNo = [];
Experiments = ["2001","2002","2003","2004","2005",...
    "2006","2007","2008","2009"];
% Experiments = ["2001","2009","2004"];
SV = struct('experimentNo',experimentNo,'Experiments',Experiments);
for exp = 1 : length(SV.Experiments)
        %Variables to update
    experimentNo = SV.Experiments(exp);
    SV.experimentNo = experimentNo; %Variables to save
    disp(experimentNo)
    if exp == 1
        StressThres = [];
        TEnerThres = [];
        SpreadEner = [];
        MPa = [];
        Table = table(MPa,StressThres,TEnerThres,SpreadEner);
        SV.Table = Table;
        SV.time2 = cell(1,length(SV.Experiments));
        SV.SpreadEner = cell(1,length(SV.Experiments));
        SV.TimeEnd = [];
        SV.Fun_Time = cell(1,length(SV.Experiments));
        SV.Fun_TensileStress = cell(1,length(SV.Experiments));
        SV.Resolution = [];
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
    SV.TimeEnd(exp) = TimeEnd;
    SV.Fun_Time{exp} = CSVDataOffs.Fun_Time;
    SV.Fun_TensileStress{exp} = CSVDataOffs.Fun_TensileStress;
    SV.Resolution(exp) = Resolution;
    
end


figure('name', 'Cumulative A-ener vs Load','Position',[60,60,1400,700])
hold on
for k = 1 : exp
%     yyaxis left
%     plot(SV.time2{k}/SV.Resolution(k), SV.SpreadEner{k});
    plot(SV.Fun_TensileStress{k}, SV.SpreadEner{k});
%     yyaxis right
%     plot(SV.Fun_Time{k}, SV.Fun_TensileStress{k})
end
title('Cumulative A-ener vs Load');
% xlim([0 max(SV.Fun_Time)]);
% xlabel('Time [s]');
xlabel('Stress [MPa]');
ylabel('Energy');
legend(num2str(round(SV.Table.MPa)));
hold off

SV.Table = sortrows(SV.Table);