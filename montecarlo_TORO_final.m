clear;
clc;

%% PARAMETERS

Eq = readmatrix("Parameters.xlsx",'Sheet','Eq(jl)','Range','B2:E21'); 
Eb = reatrix("Parameters.xlsx",'Sheet','Eb(jm)','Range','B2:D21');
EqMULT = readmatrix("Parameters.xlsx",'Sheet','EqMULT','Range','B2:E21'); 
EbMULT = readmatrix("Parameters.xlsx",'Sheet','EbMULT','Range','B2:D21');
randnum = readmatrix("Parameters.xlsx",'Sheet','RANDOM VALUES (For Monte Carlo)','Range','B1:B2000');

Eq_low = Eq .* EqMULT;
Eq_high = Eq;

Eb_low = Eb .* EbMULT;
Eb_high = Eb;

%% MONTE CARLO SIMULATION ON THE DETERMINISTIC

alpha = randnum(1:150);

Eq1 = zeros(size(Eq,1),size(Eq,2),size(alpha,1));
Eb1 = zeros(size(Eb,1),size(Eb,2),size(alpha,1));

for i = 1:size(alpha,1)
    Eq1(:,:,i) = Eq_high-alpha(i).*(Eq_high-Eq_low);
    Eq2 = squeeze(Eq1(:,:,i));

    Eb1(:,:,i) = Eb_low+alpha(i).*(Eb_high-Eb_low);
    Eb2 = squeeze(Eb1(:,:,i));

    [TOCW(i), fx(:,:,:,i), optimsol(i), fN(:,:,:,i), fO(:,:,:,i), TotCost(i), fENVI(i)] = montecarlo_function(Eq2, Eb2);
end

u = mean(TotCost);
s = std(TotCost);
n = (1.96*s/(0.05*u))^2;

cost_min = min(TotCost);
cost_max = max(TotCost);


%% TARGET SETTING

clc;

[TOCW1a, fx1a, optimsol1a, fN1a, fO1a, faa, fRa, fSa, fTa, fUa, fVa, TotCosta, fENVIa] = wwtp_TORO_function(Eq_high, Eb_high);
[TOCW1b, fx1b, optimsol1b, fN1b, fO1b, fab, fRb, fSb, fTb, fUb, fVb, TotCostb, fENVIb] = wwtp_TORO_function(Eq_high, Eb_low);
[TOCW1c, fx1c, optimsol1c, fN1c, fO1c, fac, fRc, fSc, fTc, fUc, fVc, TotCostc, fENVIc] = wwtp_TORO_function(Eq_low, Eb_high);
[TOCW1d, fx1d, optimsol1d, fN1d, fO1d, fad, fRd, fSd, fTd, fUd, fVd, TotCostd, fENVId] = wwtp_TORO_function(Eq_low, Eb_low);

bestcost = min([TotCosta TotCostb TotCostc TotCostd]);
worstcost = max([TotCosta TotCostb TotCostc TotCostd]);

A = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
budget = A*bestcost+(1-A)*worstcost;

%% TARGET-ORIENTED ROBUST OPTIMIZATION (TORO)

%Pre Allocation
Ah = ones(11,1);
Al = zeros(11,1);
theta = zeros(11,1);

for i=1:11
    while Ah(i)-Al(i)>0.05 %bisection method
        theta(i)=(Ah(i)+Al(i))/2; %midpoint computation

        Eq3(:,:,i) = Eq_high-theta(i).*(Eq_high-Eq_low);
        Eq4 = squeeze(Eq3(:,:,i));

        Eb3(:,:,i) = Eb_low+theta(i).*(Eb_high-Eb_low);
        Eb4 = squeeze(Eb3(:,:,i));

        [TOCW1(i), fx1(:,:,:,i), optimsol1(i), fN1(:,:,:,i), fO1(:,:,:,i), fa1(:,:,:,i), fR1(:,:,:,i), fS1(:,:,:,i), fT1(:,:,i), fU1(:,:,i), fV1(:,:,i), TotCost1(i), fENVI1(i)] = wwtp_TORO_function(Eq4, Eb4);

        if TotCost1(i)<=budget(i)
             Al(i)=theta(i);
        elseif TotCost1(i)>budget(i)
             Ah(i)=theta(i);
        end
    end
end

%% MONTE CARLO SIMULATION ON TORO

clc;

%MONTE CARLO ON TORO
j=rand(11,30);

for i=1:11
    for k=1:30
        Eq5(:,:,i,k) = Eq_high-j(i,k).*(Eq_high-Eq_low);
        Eq6 = squeeze(Eq5(:,:,i,k));

        Eb5(:,:,i,k) = Eb_low+j(i,k).*(Eb_high-Eb_low);
        Eb6 = squeeze(Eb5(:,:,i));

        x_MOT = squeeze(fx1(:,:,:,i));

        [TOC2(i,k), xx2(:,:,:,i,k), fval2(i,k), fN2(:,:,:,i,k), fO2(:,:,:,i,k), TotCost2(i,k), fENVI2(i,k)] = montecarlo_on_TORO_function(x_MOT, Eq6, Eb6);

        if TotCost2(i,k)<=budget(i)
            budgetsatisfied(i,k)=1;
        else
            budgetsatisfied(i,k)=0;
        end
    end
end

ave=mean(TotCost2,2);
stdev=std(TotCost2,0,2);
U=mean(budgetsatisfied,2);
sumtable = ["Theta" "Budget" "Total Cost" "Average" "Std Dev" "Total Environmental" "Probability";theta budget'/1000 TotCost1'/1000 ave/1000 stdev fENVI1'/1000 U];


%% 


% clc;
% 
% % y = (squeeze(fN(20,1,1,:)));
% % x = linspace(0,1,10);
% 
% xx = 1:10;
% yy = squeeze(fN(20,1,1,:));
% scatter(xx,yy);
% 
% hold on 
% yyy = 0.002;
% plot(xx,yyy,Color=r,LineWidth=0.5);
% hold off

%Period 1
y = (squeeze(fN(20,:,1,:)));
x = linspace(1,150,150);
scatter(x,y,'filled')
ylabel('Quality Percentage'), xlabel('Run');
% grid(ax,'on');

hold on
% yline(0.0015,'Color','r',LineWidth=2,LabelHorizontalAlignment=);
yline(0.00042587,'Color','r',LineWidth=2);
yline(0.00045714,'Color','g',LineWidth=2);
yline(0.00065871,'Color','y',LineWidth=2);
hold off

%%
%Period 1
y = (squeeze(fO(20,:,1,:)));
x = linspace(1,150,150);
scatter(x,y,'filled')
ylabel('Quality Percentage'), xlabel('Run');
% grid(ax,'on');

hold on
% yline(0.0015,'Color','r',LineWidth=2,LabelHorizontalAlignment=);
yline(0.0005,'Color','r',LineWidth=2);
yline(0.0007,'Color','g',LineWidth=2);
yline(0.0007,'Color','y',LineWidth=2);
hold off

%% 

y1 = (TotCost);
x1 = linspace(1,150,150);
scatter(x1,y1,'filled')
ylabel('Cost'), xlabel('Run');
% grid(ax,'on');

hold on
yline(mean(TotCost),'Color','y',LineWidth=2);
yline(6840.5,'Color','r',LineWidth=2);
% yline(0.0005,'Color','g',LineWidth=2);
% yline(0.005,'Color','y',LineWidth=2);
% yline(0.075,'Color','b',LineWidth=2);
hold off

%y1 = (squeeze(fO(20,:,1,:)));




% %Period 2 
% y = (squeeze(fN(20,:,2,:)));
% x = linspace(1,30,30);
% scatter(x,y,'filled')
% ylabel('Quality Percentage'), xlabel('Run');
% 
% % hold on
% % yline(0.0015,'Color','r',LineWidth=2);
% % hold off

% %Period 3
% y = (squeeze(fN(20,:,3,:)));
% x = linspace(1,30,30);
% scatter(x,y,'filled')
% ylabel('Quality Percentage'), xlabel('Run');
% 
% % hold on
% % yline(0.0015,'Color','r',LineWidth=2);
% % hold off

% %COST PLOT
% y = (optimsol);
% x = linspace(1,30,30);
% scatter(x,y,'filled')
% ylabel('COST'), xlabel('Run');
% hold on
% yline(17089.13,'Color','r',LineWidth=2);
% hold off

%% TRY-CATCH CODE (CHECKING FOR INFEASIBILITIES) FOR TORO

% for i=1:11
%     while Ah(i)-Al(i)>0.05 %bisection method
%         theta(i)=(Ah(i)+Al(i))/2; %midpoint computation
% 
%         Eq3(:,:,i) = Eq_low+theta(i).*(Eq_high-Eq_low);
%         Eq4 = squeeze(Eq3(:,:,i));
% 
%         Eb3(:,:,i) = Eb_low+theta(i).*(Eb_high-Eb_low);
%         Eb4 = squeeze(Eb3(:,:,i));
% 
%         try
%             [TOCW1(i), fx1(:,:,:,i), optimsol1(i), fN1(:,:,:,i), fO1(:,:,:,i), fa(:,:,:,i), fR(:,:,:,i), fS(:,:,:,i), fT(:,:,i), fU(:,:,i), fV(:,:,i)] = wwtp_TORO_function(Eq4, Eb4);
%         catch
%             flag = warning("Infeasible");
%         end
% 
%         if flag == "Infeasible"
%             Ah(i)=theta(i);
%         else
%             if optimsol1(i)<=budget(i)
%                  Al(i)=theta(i);
%             elseif optimsol1(i)>budget(i)
%                  Ah(i)=theta(i);
%             end
%         end
%     end
% end
