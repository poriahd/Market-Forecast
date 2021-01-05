clear
clc
close all
load('HistoricalData.mat','HisDat')

% Some definition

% The TLS folder need to be exisit in the current directory
% can be downloaded from http://www.mathworks.com/matlabcentral/fileexchange/8277, 2006.

global CrntDrctr 
CrntDrctr = pwd;

% Curve fitting data
FitCurve = {[]; 'Wind'; 'PV'; 'Wind+PV'; 'HydroPower'; 'Nuclear'; 'CHP';...
     'Thermal'; 'Net Import'; 'Supply Curve'; 'Supply Curve-(PV+Wind)';...
    'Supply Curve-(NPP+PV+Wind)'};
FitCurve(1,1:10) = {[],'All','Warm Season', 'Spring/Automn',...
    'Cold Seasons','','','','',''};

% Data selecting
TrYr = [2018 2019]; % Starting and Ending of Training data
TsYr = [2018 2019]; % Starting and Ending of Test data

% Season Definition
% Ssn = {
%             [6 7 8];
%             [3 4 5 9 10 11];
%             [12 1 2];
% };
Ssn = {[6 7 8 3 4 5 9 10 11 12 1 2];
};
%% Annual Analysis

[AnDt] = AnAnls(HisDat,TrYr);

% Area plot for annual production
Fig1 = figure('InvertHardcopy','off','Color',[1 1 1]);
Ax1 = axes('Parent',Fig1);
xlabel('Year','FontWeight','bold')
ylabel('Power (MWh)','FontWeight','bold')
xlim(Ax1,[0.95 7]);
set(Ax1,'FontSize',12,'FontWeight','bold','XTick',1:7,'XTickLabel',...
    ['2013';'2014';'2015';'2016';'2017';'2018';'2019';'2020']);
legend
hold on
Production = [AnDt.NPP(:,25) AnDt.Wind(:,25) AnDt.PV(:,25)...
    AnDt.Hyd(:,25) AnDt.CHP(:,25) AnDt.Thr(:,25) AnDt.Nimp(:,25)];    
area1 = area(Production);
set(area1(1),'DisplayName','NPP');
set(area1(2),'DisplayName','Wind');
set(area1(3),'DisplayName','PV');
set(area1(4),'DisplayName','Hyd');
set(area1(5),'DisplayName','CHP');
set(area1(6),'DisplayName','Thr');
set(area1(7),'DisplayName','Net Imp');
plot(AnDt.Cns(:,25),'k','LineWidth',2,'DisplayName','Total Demand')

% Annual Price Plot
fig2 = figure('InvertHardcopy','off','Color',[1 1 1]);
Ax2 = axes('Parent',fig2);
xlabel('Year','FontWeight','bold')
ylabel('Price (€/MWh)','FontWeight','bold')
xlim(Ax2,[1 7]);
ylim(Ax2,[20 50]);
set(Ax2,'FontSize',10,'FontWeight','bold','XTick',1:7,'XTickLabel',...
    ['2013';'2014';'2015';'2016';'2017';'2018';'2019';'2020']);
hold on
plot(AnDt.DyAhPr(:,25),'k','LineWidth',2) 

% Annual Corelation analysis 
AnPrc = VeNr(AnDt.DyAhPr(:,25)); % Normalise the values
AnPrdNr = [VeNr(AnDt.NPP(:,25)) VeNr(AnDt.Wind(:,25)) VeNr(AnDt.PV(:,25))...
    VeNr(AnDt.Hyd(:,25)) VeNr(AnDt.CHP(:,25)) VeNr(AnDt.Thr(:,25)) ...
    VeNr(AnDt.Nimp(:,25)) VeNr(AnDt.Tprd(:,25)) VeNr(AnDt.Cns(:,25))];

[RhoAn,PvalAn] = corr(AnPrc,AnPrdNr);

% Corelation analysis price with different productions
[SlDt,~] = DataSelector(HisDat,2013,2019,1:12);

Prc = VeNr(SlDt.DyAhPr); %circshift(Prc,1)
PrdNr = [VeNr(SlDt.NPP) VeNr(SlDt.Wind) VeNr(SlDt.PV) ...
    VeNr(SlDt.Hyd) VeNr(SlDt.CHP) VeNr(SlDt.Thr) ...
    VeNr(SlDt.Nimp) VeNr(SlDt.Tprd) VeNr(SlDt.Cns)];
    
[Rho,Pval] = corr(Prc,PrdNr);

% Corelation analysis price with different combination of demand & RES/NPP
[SlDt,~] = DataSelector(HisDat,2013,2019,1:12);
Prc = VeNr(SlDt.DyAhPr); 
SupNr = [VeNr(SlDt.Wind+SlDt.PV) ...
     VeNr(SlDt.Cns) VeNr(SlDt.Cns-SlDt.Wind-SlDt.PV)...
     VeNr(SlDt.Cns-SlDt.Wind-SlDt.PV-SlDt.NPP)...
     ];
[RhoS,PvalS] = corr(Prc,SupNr);


%% Curve price/Production type

[SlDt, ~] = DataSelector(HisDat,2017,2019,1:12);

fig3 = figure('InvertHardcopy','off','Color',[1 1 1]);
PltDv = [2 4]; % subplot devision m*n

% Wind Curve
TmpPlt = {'Wind' PltDv 1};
[Ffit] = PriceElacity(SlDt.Wind,SlDt.DyAhPr,TmpPlt);
FitCurve{2,2} = [Ffit.p1 Ffit.p2];

% PV Curve
TmpPlt = {'PV' PltDv 2};
[Ffit] = PriceElacity(SlDt.PV,SlDt.DyAhPr,TmpPlt);
FitCurve{3,2} = [Ffit.p1 Ffit.p2];

% Wind+PV Curve
TmpPlt = {'Wind + PV' PltDv 3};
[Ffit] = PriceElacity(SlDt.Wind + SlDt.PV,SlDt.DyAhPr,TmpPlt);
FitCurve{4,2} = [Ffit.p1 Ffit.p2];

% Hydro Curve
TmpPlt = {'Hydro' PltDv 4};
[Ffit] = PriceElacity(SlDt.Hyd,SlDt.DyAhPr,TmpPlt);
FitCurve{5,2} = [Ffit.p1 Ffit.p2];

% NPP Curve
TmpPlt = {'NPP' PltDv 5};
[Ffit] = PriceElacity(SlDt.NPP,SlDt.DyAhPr,TmpPlt);
FitCurve{6,2} = [Ffit.p1 Ffit.p2];
ylabel({'Price (€/MWh)';' '},'FontWeight','bold','FontSize',12) 

% CHP Curve
TmpPlt = {'CHP' PltDv 6};
[Ffit] = PriceElacity(SlDt.CHP,SlDt.DyAhPr,TmpPlt);
FitCurve{7,2} = [Ffit.p1 Ffit.p2];

% Thermal Curve
TmpPlt = {'Other Thermal' PltDv 7};
[Ffit] = PriceElacity(SlDt.Thr,SlDt.DyAhPr,TmpPlt);
FitCurve{8,2} = [Ffit.p1 Ffit.p2];

% Net Import Curve
TmpPlt = {'Net Import' PltDv 8};
[Ffit] = PriceElacity(SlDt.Nimp,SlDt.DyAhPr,TmpPlt);
FitCurve{9,2} = [Ffit.p1 Ffit.p2];
xlabel({' ';'Production (MW)'},'FontWeight','bold','FontSize',12)


%% Calculate Supply curve

fig4 = figure('InvertHardcopy','off','Color',[1 1 1]);
PltDv = [size(Ssn,1) 3]; % subplot devision m*n

for SsnCnt = 1:size(Ssn,1)
    
    [TrDt, ~] = DataSelector(HisDat,TrYr(1,1),TrYr(1,2),Ssn{SsnCnt});

    % Total Supply Curve (Method 1)
    TmpPlt = {'Method 1' PltDv 3*(SsnCnt-1)+1};
    [Ffit] = SuplyCurve(TrDt.DyAhPr,TrDt.Cns,TmpPlt);
    FitCurve{10,SsnCnt+2} = Ffit;
    
    % Supply Curve - Wind - PV (Method 2)
    TmpPlt = {'Method 2' PltDv 3*(SsnCnt-1)+2};
    [Ffit] = SuplyCurve(TrDt.DyAhPr,TrDt.Cns - ...
        TrDt.Wind - TrDt.PV,TmpPlt);
    FitCurve{11,SsnCnt+2} = Ffit;
    
    % Supply Curve - Wind - PV - NPP (Method 3)
    TmpPlt = {'Method 3' PltDv 3*(SsnCnt-1)+3};
    [Ffit] = SuplyCurve(TrDt.DyAhPr,TrDt.Cns - ...
        TrDt.Wind - TrDt.PV - TrDt.NPP,TmpPlt);
    FitCurve{12,SsnCnt+2} = Ffit; 
end

save('FitREs.mat','FitCurve')

%% Test the Estimation

EsTsPrC = cell(size(Ssn,1),3);
TsEr = zeros(size(Ssn,1),3);
PltClr = {'k';'g';'r';'m';'c';'b';'y';};
MthdNm = {'Method 1';'Method 2';'Method 3'};

fig5 = figure('InvertHardcopy','off','Color',[1 1 1]);
PltDv = [1 size(Ssn,1)]; % subplot devision m*n

for SsnCnt = 1:size(Ssn,1)
    
    [TsDt, ~] = DataSelector(HisDat,TsYr(1,1),TsYr(1,2),Ssn{SsnCnt});
    
    % calculate average profile for estimation
    AvTsDt.DyAhPr(SsnCnt,:) = mean(TsDt.DyAhPr);        
    AvTsDt.SsnNm(SsnCnt,1) = SsnCnt;
    AvTsDt.Tprd(SsnCnt,:) = mean(TsDt.Tprd);
    AvTsDt.Hyd(SsnCnt,:) = mean(TsDt.Hyd);
    AvTsDt.NPP(SsnCnt,:) = mean(TsDt.NPP);
    AvTsDt.Wind(SsnCnt,:) = mean(TsDt.Wind);
    AvTsDt.CHP(SsnCnt,:) = mean(TsDt.CHP);
    AvTsDt.Thr(SsnCnt,:) = mean(TsDt.Thr);
    AvTsDt.PV(SsnCnt,:) = mean(TsDt.PV);
    AvTsDt.Nimp(SsnCnt,:) = mean(TsDt.Nimp);
    AvTsDt.Cns(SsnCnt,:) = mean(TsDt.Cns);
    
    Ax5 = subplot(PltDv(1),PltDv(2),SsnCnt);
    plot(AvTsDt.DyAhPr(SsnCnt,:),PltClr{1},'LineWidth',1,...
                'DisplayName','Average Price')        
    title(['Season ' num2str(SsnCnt)]);
    xlabel('Time (h)','FontWeight','bold')
    ylabel('Price (€/MWh)','FontWeight','bold')
    xlim(Ax5,[0 24]);
    set(Ax5,'FontSize',10,'FontWeight','bold')
    ylim([20 60])
    legend
    hold on
    
    % Estimate the average Price in Season:
    % Total Supply Curve (Method 1)
    EsTsPrC{SsnCnt,1} = polyval(FitCurve{10,SsnCnt+2},...
        AvTsDt.Cns(SsnCnt,:));

    % Supply Curve - Wind - PV (Method 2)
    EsTsPrC{SsnCnt,2} = polyval(FitCurve{11,SsnCnt+2},...
        AvTsDt.Cns(SsnCnt,:) - AvTsDt.Wind(SsnCnt,:) - ...
        AvTsDt.PV(SsnCnt,:));
    
    % Supply Curve - Wind - PV - NPP (Method 3)
    EsTsPrC{SsnCnt,3} = polyval(FitCurve{12,SsnCnt+2},...
        AvTsDt.Cns(SsnCnt,:) - AvTsDt.Wind(SsnCnt,:) - ...
        AvTsDt.PV(SsnCnt,:) - AvTsDt.NPP(SsnCnt,:));
  
    % Error calculation        
    for MtdCnt = 1:3
        TsEr(SsnCnt,MtdCnt) = mean(abs(EsTsPrC{SsnCnt,MtdCnt} - ...
            AvTsDt.DyAhPr(SsnCnt,:))./abs(AvTsDt.DyAhPr(SsnCnt,:)))*100;
      
        plot(EsTsPrC{SsnCnt,MtdCnt},PltClr{MtdCnt+1},...
                'LineWidth',1,'DisplayName',MthdNm{MtdCnt})
    end    
end

%% Estimate 2030
FctNPP = 1.57;
FctWind = 2.7;
FctPV = 10;
FctCns1 = 1.1;
FctCns2 = 1.2;

fig6 = figure('InvertHardcopy','off','Color',[1 1 1]);
PltDv = [2 size(Ssn,1)]; % subplot devision m*n

for SsnCnt = 1:size(Ssn,1)
    % Base case :2018
    [Dt2018, ~] = DataSelector(HisDat,2018,2018,Ssn{SsnCnt});
    
    % calculate average profile for estimation
    Dt2030.SsnNm(SsnCnt,1) = SsnCnt;
    
    Dt2030.NPP(SsnCnt,:) = FctNPP*mean(Dt2018.NPP);
    Dt2030.Wind(SsnCnt,:) = FctWind*mean(Dt2018.Wind);
    Dt2030.PV(SsnCnt,:) = FctPV*mean(Dt2018.PV);

% Scenario #1
    Dt2030.Cns(SsnCnt,:) = FctCns1*mean(Dt2018.Cns)
    
    Ax6 = subplot(PltDv(1),PltDv(2),SsnCnt);
    title(['Season ' num2str(SsnCnt)]);
    ylabel('Price (€/MWh)','FontWeight','bold')
    xlim(Ax6,[0 24]);
    set(Ax6,'FontSize',10,'FontWeight','bold')
    legend
    hold on
    
    % Estimate the average Price in Season:
    % Total Supply Curve (Method 1)
    EsPrC2030S1{SsnCnt,1} = polyval(FitCurve{10,SsnCnt+2},...
        Dt2030.Cns(SsnCnt,:));    
    plot(EsPrC2030S1{SsnCnt,1},PltClr{1},'LineWidth',1,...
        'DisplayName','Method 1')
            
    % Supply Curve - Wind - PV (Method 2)
    EsPrC2030S1{SsnCnt,2} = polyval(FitCurve{11,SsnCnt+2},...
        Dt2030.Cns(SsnCnt,:) - Dt2030.Wind(SsnCnt,:) - ...
        Dt2030.PV(SsnCnt,:));
    plot(EsPrC2030S1{SsnCnt,2},PltClr{2},'LineWidth',1,...
        'DisplayName','Method 2')
    
    % Supply Curve - Wind - PV - NPP (Method 3)
    EsPrC2030S1{SsnCnt,3} = polyval(FitCurve{12,SsnCnt+2},...
        Dt2030.Cns(SsnCnt,:) - Dt2030.Wind(SsnCnt,:) - ...
        Dt2030.PV(SsnCnt,:) - Dt2030.NPP(SsnCnt,:));
    plot(EsPrC2030S1{SsnCnt,3},PltClr{3},'LineWidth',1,...
        'DisplayName','Method 3')

    Dt2030.Cns(SsnCnt,:)
% Scenario #2
    Dt2030.Cns(SsnCnt,:) = FctCns2*mean(Dt2018.Cns);
    
    Ax6 = subplot(PltDv(1),PltDv(2),size(Ssn,1)+SsnCnt);
    title(['Season ' num2str(SsnCnt)]);
    xlabel('Time (h)','FontWeight','bold')
    xlim(Ax6,[0 24]);
    set(Ax6,'FontSize',10,'FontWeight','bold')
    legend
    hold on
    
    % Estimate the average Price in Season:
    % Total Supply Curve (Method 1)
    EsPrC2030S2{SsnCnt,1} = polyval(FitCurve{10,SsnCnt+2},...
        Dt2030.Cns(SsnCnt,:));    
    plot(EsPrC2030S2{SsnCnt,1},PltClr{1},'LineWidth',1,...
        'DisplayName','Method 1')
            
    % Supply Curve - Wind - PV (Method 2)
    EsPrC2030S2{SsnCnt,2} = polyval(FitCurve{11,SsnCnt+2},...
        Dt2030.Cns(SsnCnt,:) - Dt2030.Wind(SsnCnt,:) - ...
        Dt2030.PV(SsnCnt,:));
    plot(EsPrC2030S2{SsnCnt,2},PltClr{2},'LineWidth',1,...
        'DisplayName','Method 2')
    
    % Supply Curve - Wind - PV - NPP (Method 3)
    EsPrC2030S2{SsnCnt,3} = polyval(FitCurve{12,SsnCnt+2},...
        Dt2030.Cns(SsnCnt,:) - Dt2030.Wind(SsnCnt,:) - ...
        Dt2030.PV(SsnCnt,:) - Dt2030.NPP(SsnCnt,:));
    plot(EsPrC2030S2{SsnCnt,3},PltClr{3},'LineWidth',1,...
        'DisplayName','Method 3')
    
end
%% Functions: 

function [AnDt] = AnAnls(HisDat,TrYr)
% Make Annual Profiles

    for AnCnt = 1:TrYr(2)-TrYr(1)+1
        Year = AnCnt-1+TrYr(1);
        AnDt.Date(AnCnt,:) = Year;
        [SlDt,~] = DataSelector(HisDat,Year,Year,1:12);
        
        AnDt.DyAhPr(AnCnt,1:24) = mean(SlDt.DyAhPr);
        AnDt.Tprd(AnCnt,1:24) = sum(SlDt.Tprd);
        AnDt.Hyd(AnCnt,1:24) = sum(SlDt.Hyd);
        AnDt.NPP(AnCnt,1:24) = sum(SlDt.NPP);
        AnDt.Wind(AnCnt,1:24) = sum(SlDt.Wind);
        AnDt.CHP(AnCnt,1:24) = sum(SlDt.CHP);
        AnDt.Thr(AnCnt,1:24) = sum(SlDt.Thr);
        AnDt.PV(AnCnt,1:24) = sum(SlDt.PV);
        AnDt.Nimp(AnCnt,1:24) = sum(SlDt.Nimp);
        AnDt.Tprd(AnCnt,1:24) = sum(SlDt.Tprd);
        AnDt.Cns(AnCnt,1:24) = sum(SlDt.Cns);
    
        AnDt.DyAhPr(AnCnt,25) = mean(AnDt.DyAhPr(AnCnt,1:24));
        AnDt.Tprd(AnCnt,25) = sum(AnDt.Tprd(AnCnt,1:24));
        AnDt.Hyd(AnCnt,25) = sum(AnDt.Hyd(AnCnt,1:24));
        AnDt.NPP(AnCnt,25) = sum(AnDt.NPP(AnCnt,1:24));
        AnDt.Wind(AnCnt,25) = sum(AnDt.Wind(AnCnt,1:24));
        AnDt.CHP(AnCnt,25) = sum(AnDt.CHP(AnCnt,1:24));
        AnDt.Thr(AnCnt,25) = sum(AnDt.Thr(AnCnt,1:24));
        AnDt.PV(AnCnt,25) = sum(AnDt.PV(AnCnt,1:24));
        AnDt.Nimp(AnCnt,25) = sum(AnDt.Nimp(AnCnt,1:24));
        AnDt.Tprd(AnCnt,25) = sum(AnDt.Tprd(AnCnt,1:24));
        AnDt.Cns(AnCnt,25) = sum(AnDt.Cns(AnCnt,1:24));
    end
end

% Selecting data
function [SlDt, Ind] = DataSelector(HisDat,StYr,EnYr,MntRange)
    switch StYr % 2010:1 2011:366 2012:731
        case 2013
            St = 1097;
        case 2014
            St = 1462;
        case 2015
            St = 1827;
        case 2016
            St = 2192;
        case 2017
            St = 2558;
        case 2018
            St = 2923;
        case 2019
            St = 3288;
    end
    switch EnYr
        case 2013
            En = 1461;
        case 2014
            En = 1826;
        case 2015
            En = 2191;
        case 2016
            En = 2557;
        case 2017
            En = 2922;
        case 2018
            En = 3287;
        case 2019
            En = 3652;   
    end
    
    Ind = false(size(HisDat.Date,1),1);
    for Cnt = 1: size(HisDat.Date,1)
        if Cnt >= St && Cnt <= En
            tempmonth = month(datetime(HisDat.Date(Cnt,1)));
            if ismember(tempmonth,MntRange)
                Ind(Cnt,1) = true;
            end
        end
    end
    SlDt.Date = HisDat.Date(Ind);
    SlDt.DyAhPr = HisDat.DyAhPr(Ind,:);
    SlDt.Tprd = HisDat.Tprd(Ind,:);
    SlDt.Hyd = HisDat.Hyd(Ind,:);
    SlDt.NPP = HisDat.NPP(Ind,:);
    SlDt.Wind = HisDat.Wind(Ind,:);
    SlDt.CHP = HisDat.CHP(Ind,:);
    SlDt.Thr = HisDat.Thr(Ind,:);
    SlDt.PV = HisDat.PV(Ind,:);
    SlDt.Nimp = HisDat.Nimp(Ind,:);
    SlDt.Cns = HisDat.Cns(Ind,:);  
    
end

function [y] = VeNr(x)
    % Reshape to a vector and Normalise X
    temp = reshape(x,[],1);
    y = temp./max(temp);    
end

function [Ffit] = PriceElacity(x,y,TmpPlt)
    % Create data for elacity curve
    
    X = reshape(x,[],1);
    Y = reshape(y,[],1);        
    [~,Iy] = rmoutliers(Y,'median'); % Clean outlier
    X = X(~Iy);
    Y = Y(~Iy);    
    [WX,WY,W] = GroupingData(X,Y);
    Ffit = fit(X,Y,'poly1');
        
    Ax3 = subplot(TmpPlt{2}(1),TmpPlt{2}(2),TmpPlt{3});
    set(Ax3,'FontSize',10,'YTick',[0 40 80],'FontWeight','bold')
    title(TmpPlt{1},'FontSize',10);
    grid on
    ylim([0 80])
    hold on

    scatter(WX,WY,3*W,'Marker','.','MarkerEdgeColor',[1 0 0])

    Xf = linspace(0,max(WX));
    Yf = Ffit(Xf);
    plot(Xf,Yf,'k','LineWidth',1.5) 
    
end

function [X,Y,W] = GroupingData(x,y)
    X = [];
    Y = [];
    W = [];
    [N,c] = hist3([x,y],[1000,500]);
    for Cnt1 = 1:size(N,1)
        for Cnt2 = 1:size(N,2)
            if N(Cnt1,Cnt2) ~= 0
               X = [X; c{1,1}(1,Cnt1)];
               Y = [Y; c{1,2}(1,Cnt2)];
               W = [W; N(Cnt1,Cnt2)];
            end
        end
    end
end

function [Ffit] = SuplyCurve(x,y,TmpPlt)
    % Create data for elacity curve
    global CrntDrctr

    X = reshape(x,[],1);
    Y = reshape(y,[],1);        
     
    % Clean outlier
    [~,Ix] = rmoutliers(X,'median');
    [~,Iy] = rmoutliers(Y,'median');
    X = X(~Ix & ~Iy);
    Y = Y(~Ix & ~Iy);    
    [WX,WY,W] = GroupingData(X,Y);

    cd([CrntDrctr '/TLS'])
    [~, Ffit] = fit_2D_data(Y, X, 'no');
    cd(CrntDrctr)
    
    Ax4 = subplot(TmpPlt{2}(1),TmpPlt{2}(2),TmpPlt{3});
    set(Ax4,'FontSize',10,'FontWeight','bold')
    Ax4.YAxis.Exponent = 3;
    xlim(Ax4,[10 70]);
    ylim(Ax4,[1000*floor(min(Y)/1000) 1000*ceil(max(Y)/1000)]);
    grid on
    hold on

    scatter(WX,WY,3*W,'Marker','.','MarkerEdgeColor',[1 0 0])
        
    Yf = linspace(min(Y),max(Y));
    Xf = polyval(Ffit,Yf);
    hold on
    plot(Xf,Yf,'k','LineWidth',1)

    if TmpPlt{3}<=3
        title({TmpPlt{1};' '},'FontSize',10);
    end
    if TmpPlt{3} == TmpPlt{2}(1)*TmpPlt{2}(2)-1
        xlabel({' ';'Price (€/MWh)'},'FontWeight','bold')
    end
    if mod(TmpPlt{3},TmpPlt{2}(1)) == 1
        if ceil(TmpPlt{3}/TmpPlt{2}(1)) == ceil(TmpPlt{2}(1)/2)
            ylabel({['Season ' num2str(floor(TmpPlt{3}/TmpPlt{2}(1))+1)]...
                ;' '; 'Power (MW)'});
        else
            ylabel({['Season ' ...
                num2str(floor(TmpPlt{3}/TmpPlt{2}(1))+1)];' ';' '});
        end            
    end            
end    


