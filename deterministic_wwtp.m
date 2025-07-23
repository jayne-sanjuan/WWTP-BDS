clear
clc

%INDECES
j = 20; %nodes
l = 4; %quality parameter
m = 3; %byproduct
n = 5; %waste processing
o = 2; %bioproduct
p = 3; %landfill
q = 3; %agri soil
r = 2; %turbine generator
s = 4; %storage
t = 3; %period
u = 4; %"to" clusters
v = 5; %"from" clusters

%PARAMETERS NETWORK
f = readmatrix("Parameters.xlsx",'Sheet','f(it)','Range','B2:D4'); %water input
cm = readmatrix("Parameters.xlsx",'Sheet','cm(jk)','Range','B2:B13'); %capacity for machine k
cs = readmatrix("Parameters.xlsx",'Sheet','cs(s)','Range','B1:B4'); %capacity for storage s
dw = readmatrix("Parameters.xlsx",'Sheet','dw(t)','Range','B2:B4'); %demand of water
M = [100100; 112000; 400000]; %M = max(cm,[],"all");

%PARAMETERS EFFECTIVENESS
L(:,:,1) = readmatrix("Parameters.xlsx",'Sheet','L(ilt)','Range','B2:E4');
L(:,:,2) = readmatrix("Parameters.xlsx",'Sheet','L(ilt)','Range','H2:K4'); 
L(:,:,3) = readmatrix("Parameters.xlsx",'Sheet','L(ilt)','Range','N2:Q4'); 
K(:,:,1) = readmatrix("Parameters.xlsx",'Sheet','K(imt)','Range','B2:D4'); 
K(:,:,2) = readmatrix("Parameters.xlsx",'Sheet','K(imt)','Range','G2:I4'); 
K(:,:,3) = readmatrix("Parameters.xlsx",'Sheet','K(imt)','Range','L2:N4'); 
Eq = readmatrix("Parameters.xlsx",'Sheet','Eq(jl)','Range','B2:E21'); 
Eb = readmatrix("Parameters.xlsx",'Sheet','Eb(jm)','Range','B2:D21'); 
P = readmatrix("Parameters.xlsx",'Sheet','P(l)','Range','B2:B5'); 
W = readmatrix("Parameters.xlsx",'Sheet','W(m)','Range','B2:B4'); 
domJ = [zeros(4,3); 1 0 0; 1 0 0; 1 0 0; zeros(1,3); 0 1 0; 0 1 0; 0 1 0; zeros(1,3); 0 0 1; 0 0 1; 0 0 1; zeros(5,3)];
domR = [0 1 0 0 1; 1 1 0 0 0; 1 1 0 0 0];

%WWTP COSTS
OCw = readmatrix("Parameters.xlsx",'Sheet','OCw(jj)','Range','B2:U21');
TCw = readmatrix("Parameters.xlsx",'Sheet','TCw(jj)','Range','B2:U21');
SC = readmatrix("Parameters.xlsx",'Sheet','SC(jj)','Range','B2:U21');
HC = readmatrix("Parameters.xlsx",'Sheet','HC(j)','Range','B2:B21');
SPw = readmatrix("Parameters.xlsx",'Sheet','SPw(jj)','Range','B2:U21');

%WWTP ENV
PEw = readmatrix("Parameters.xlsx",'Sheet','PEw(jj)','Range','B2:U21');
TEw = readmatrix("Parameters.xlsx",'Sheet','TEw(jj)','Range','B2:U21');

%PARAMETERS BDS
dl = readmatrix("Parameters.xlsx",'Sheet','dl(pt)','Range','B2:D4');
dc = readmatrix("Parameters.xlsx",'Sheet','dc(qt)','Range','B2:D4');
du = readmatrix("Parameters.xlsx",'Sheet','du(t)','Range','B2:B4');
cp = readmatrix("Parameters.xlsx",'Sheet','cp(n)','Range','B2:B6'); 
Q = readmatrix("Parameters.xlsx",'Sheet','Q(o)','Range','B2:B3');
C = readmatrix("Parameters.xlsx",'Sheet','C','Range','B2:B2');
D = readmatrix("Parameters.xlsx",'Sheet','D(r)','Range','B2:B3');

%BDS COSTS
OCb = readmatrix("Parameters.xlsx",'Sheet','OCb(mnt)','Range','B2:F4');
OCn = readmatrix("Parameters.xlsx",'Sheet','OCn(ont)','Range','B2:F3');
OCt = readmatrix("Parameters.xlsx",'Sheet','OCt(rt)','Range','B2:D3');
TCb = readmatrix("Parameters.xlsx",'Sheet','TCb(mnt)','Range','B2:F4');
TCn = readmatrix("Parameters.xlsx",'Sheet','TCn(ont)','Range','B2:F3');
TCl = readmatrix("Parameters.xlsx",'Sheet','TCl(pt)','Range','B2:D4');
TCc = readmatrix("Parameters.xlsx",'Sheet','TCc(qt)','Range','B2:D4');
DCb = readmatrix("Parameters.xlsx",'Sheet','DCb(mnt)','Range','B2:F4');
DCl = readmatrix("Parameters.xlsx",'Sheet','DCl(pt)','Range','B2:D4');
SPc = readmatrix("Parameters.xlsx",'Sheet','SPc(qt)','Range','B2:D4');

%BDS ENV
PEb = readmatrix("Parameters.xlsx",'Sheet','PEb(mn)','Range','B2:F4');
DEb = readmatrix("Parameters.xlsx",'Sheet','DEb(mn)','Range','B2:F4');
PEn = readmatrix("Parameters.xlsx",'Sheet','PEn(on)','Range','B2:F3');
DEn = readmatrix("Parameters.xlsx",'Sheet','DEn(on)','Range','B2:F3');
PEt = readmatrix("Parameters.xlsx",'Sheet','PEt(r)','Range','B1:B2');
DEl = readmatrix("Parameters.xlsx",'Sheet','DEl(p)','Range','B1:B3');
PEc = readmatrix("Parameters.xlsx",'Sheet','PEc(q)','Range','B1:B3');
GHG = readmatrix("Parameters.xlsx",'Sheet','GHG(on)','Range','B2:F3');

%DECISION VARIABLES WWTP
a = optimvar('a', j, j, t, 'Type', 'integer','LowerBound', 0, 'UpperBound', 1);
b = optimvar('b', j, u, t, 'Type', 'integer','LowerBound', 0, 'UpperBound', 1);
e = optimvar('e', v, j, t, 'Type', 'integer','LowerBound', 0, 'UpperBound', 1);
x = optimvar('x', j, j, t, 'LowerBound', 0, 'UpperBound', inf);

%DECISION VARIABLES BDS
R = optimvar('R', m, n, t, 'LowerBound', 0, 'UpperBound', inf);
S = optimvar('S', o, n, t, 'LowerBound', 0, 'UpperBound', inf);
T = optimvar('T', p, t, 'LowerBound', 0, 'UpperBound', inf);
U = optimvar('U', q, t, 'LowerBound', 0, 'UpperBound', inf);
V = optimvar('V', r, t, 'LowerBound', 0, 'UpperBound', inf);

%SYSTEM VARIABLES WWTP
Ib = optimvar('Ib', j, t,'LowerBound', 0, 'UpperBound', inf); %beginning inventory
Ie = optimvar('Ie', j, t, 'LowerBound', 0, 'UpperBound', inf); %ending inventory
F = optimvar('F', j, l, t, 'LowerBound', 0, 'UpperBound', inf); %volume of remaining quality parameter in the end of treatment
G = optimvar('G', j, m, t, 'LowerBound', 0, 'UpperBound', inf); %volume of remaining byproducts in the end of treatment
H = optimvar('H', j, t, 'LowerBound', 0, 'UpperBound', inf); %volume of remaining clean water
J = optimvar('J', j, m, t, 'LowerBound', 0, 'UpperBound', inf); %volume of byproducts FOR BDS
N = optimvar('N', j, l, t, 'LowerBound', 0, 'UpperBound', inf); %remaining percent quality parameter in the end of treatment
O = optimvar('O', j, m, t, 'LowerBound', 0, 'UpperBound', inf); %remaining percent byproducts in the end of treatment
JJ = optimvar('JJ', m, t, 'LowerBound', 0, 'UpperBound', inf); %%dump J
RR = optimvar('RR', m, n, t, 'LowerBound', 0, 'UpperBound', inf);


%SYSTEM VARIABLES BDS
B = optimvar('B', o, t, 'Type', 'continuous');
A = optimvar('A', t, 'Type', 'continuous');
udl = optimvar('udl', p, t, 'LowerBound', 0, 'UpperBound', inf);
udc = optimvar('udc', q, t, 'LowerBound', 0, 'UpperBound', inf);
udu = optimvar('udu', t, 'LowerBound', 0, 'UpperBound', inf);

%OBJ VARIABLE
ECON = optimvar('ECON', 'Type', 'continuous');
ENVI = optimvar('ENVI', 'Type', 'continuous');

%RANDOM VALUES
randnum = readmatrix("Parameters.xlsx",'Sheet','RANDOM VALUES (For Monte Carlo)','Range','B1:B2000');

NETWORK = optimproblem;

NETWORK.Constraints.neta = e(1,4,:) == 1;
NETWORK.Constraints.netb = e(1,4,:) >= b(4,1,:);
NETWORK.Constraints.netc = b(4,1,:) >= e(2,8,:);
NETWORK.Constraints.netd = e(2,8,:) >= b(8,2,:) + b(8,3,:) + b(8,4,:);
NETWORK.Constraints.nete = b(8,2,:) >= e(3,12,:);
NETWORK.Constraints.netf = e(3,12,:) >= b(12,3,:) + b(12,4,:);
NETWORK.Constraints.netg = b(8,3,:) + b(12,3,:) >= e(4,16,:);
NETWORK.Constraints.neth = e(4,16,:) >= b(16,4,:);
NETWORK.Constraints.neti = b(8,4,:) + b(12,4,:) + b(16,4,:) >= e(5,20,:);
NETWORK.Constraints.netj = e(5,20,:) >= 1;

%clustering (to treatments) constraints 
NETWORK.Constraints.tocluster1 = b(4,1,:) >= a(4,5,:);
NETWORK.Constraints.tocluster2 = b(4,1,:) >= a(4,6,:);
NETWORK.Constraints.tocluster3 = b(4,1,:) >= a(4,7,:);
NETWORK.Constraints.tocluster4 = b(8,2,:) >= a(8,9,:);
NETWORK.Constraints.tocluster5 = b(8,2,:) >= a(8,10,:);
NETWORK.Constraints.tocluster6 = b(8,2,:) >= a(8,11,:);
NETWORK.Constraints.tocluster7 = b(8,3,:) >= a(8,13,:);
NETWORK.Constraints.tocluster8 = b(8,3,:) >= a(8,14,:);
NETWORK.Constraints.tocluster9 = b(8,3,:) >= a(8,15,:);
NETWORK.Constraints.tocluster10 = b(8,4,:) >= a(8,17,:);
NETWORK.Constraints.tocluster11 = b(8,4,:) >= a(8,18,:);
NETWORK.Constraints.tocluster12 = b(8,4,:) >= a(8,19,:);
NETWORK.Constraints.tocluster13 = b(12,3,:) >= a(12,13,:);
NETWORK.Constraints.tocluster14 = b(12,3,:) >= a(12,14,:);
NETWORK.Constraints.tocluster15 = b(12,3,:) >= a(12,15,:);
NETWORK.Constraints.tocluster16 = b(12,4,:) >= a(12,17,:);
NETWORK.Constraints.tocluster17 = b(12,4,:) >= a(12,18,:);
NETWORK.Constraints.tocluster18 = b(12,4,:) >= a(12,19,:);
NETWORK.Constraints.tocluster19 = b(16,4,:) >= a(16,17,:);
NETWORK.Constraints.tocluster20 = b(16,4,:) >= a(16,18,:);
NETWORK.Constraints.tocluster21 = b(16,4,:) >= a(16,19,:);
NETWORK.Constraints.tocluster22 = sum(b(4,:,:),2) == 1;
NETWORK.Constraints.tocluster23 = sum(b(8,:,:),2) == 1;
NETWORK.Constraints.tocluster24 = sum(b(12,:,:),2) == 1;
NETWORK.Constraints.tocluster25 = sum(b(16,:,:),2) == 1;

%clustering constraints (to storage)
NETWORK.Constraints.fromcluster1 = e(1,4,:) >= a(1,4,:);
NETWORK.Constraints.fromcluster2 = e(1,4,:) >= a(2,4,:);
NETWORK.Constraints.fromcluster3 = e(1,4,:) >= a(3,4,:);
NETWORK.Constraints.fromcluster4 = e(2,8,:) >= a(5,8,:);
NETWORK.Constraints.fromcluster5 = e(2,8,:) >= a(6,8,:);
NETWORK.Constraints.fromcluster6 = e(2,8,:) >= a(7,8,:);
NETWORK.Constraints.fromcluster7 = e(3,12,:) >= a(9,12,:);
NETWORK.Constraints.fromcluster8 = e(3,12,:) >= a(10,12,:);
NETWORK.Constraints.fromcluster9 = e(3,12,:) >= a(11,12,:);
NETWORK.Constraints.fromcluster10 = e(4,16,:) >= a(13,16,:);
NETWORK.Constraints.fromcluster11 = e(4,16,:) >= a(14,16,:);
NETWORK.Constraints.fromcluster12 = e(4,16,:) >= a(15,16,:);
NETWORK.Constraints.fromcluster13 = e(5,20,:) >= a(17,20,:);
NETWORK.Constraints.fromcluster14 = e(5,20,:) >= a(18,20,:);
NETWORK.Constraints.fromcluster15 = e(5,20,:) >= a(19,20,:);
NETWORK.Constraints.fromcluster16 = sum(e(:,4,:),1) == 1;
NETWORK.Constraints.fromcluster17 = sum(e(:,8,:),1) == 1;
NETWORK.Constraints.fromcluster18 = sum(e(:,12,:),1) == 1;
NETWORK.Constraints.fromcluster19 = sum(e(:,16,:),1) == 1;
NETWORK.Constraints.fromcluster20 = sum(e(:,20,:),1) == 1;

%NEW
NETWORK.Constraints.binactive1 = x(:,:,1) <= M(1)*a(:,:,1);
NETWORK.Constraints.binactive2 = x(:,:,2) <= M(2)*a(:,:,2);
NETWORK.Constraints.binactive3 = x(:,:,3) <= M(3)*a(:,:,3);

%input
NETWORK.Constraints.input1 = f(1,:)' == squeeze(x(1,4,:));
NETWORK.Constraints.input2 = f(2,:)' == squeeze(x(2,4,:));
NETWORK.Constraints.input3 = f(3,:)' == squeeze(x(3,4,:));

%demand 
NETWORK.Constraints.demand1 = squeeze(x(17,20,:) + x(18,20,:) + x(19,20,:)) >= dw;

%flow
NETWORK.Constraints.flow1 = x(1,4,:) + x(2,4,:) + x(3,4,:) == x(4,5,:) + x(4,6,:) + x(4,7,:);
NETWORK.Constraints.flow2 = x(4,5,:) >= x(5,8,:);
NETWORK.Constraints.flow3 = x(4,6,:) >= x(6,8,:);
NETWORK.Constraints.flow4 = x(4,7,:) >= x(7,8,:);
NETWORK.Constraints.flow5 = x(5,8,:) + x(6,8,:) + x(7,8,:) == x(8,9,:) + x(8,10,:) + x(8,11,:) + x(8,13,:) + x(8,14,:) + x(8,15,:) + x(8,17,:) + x(8,18,:) + x(8,19,:);
NETWORK.Constraints.flow6 = x(8,9,:) >= x(9,12,:);
NETWORK.Constraints.flow7 = x(8,10,:) >= x(10,12,:);
NETWORK.Constraints.flow8 = x(8,11,:) >= x(11,12,:);
NETWORK.Constraints.flow9 = x(9,12,:) + x(10,12,:) + x(11,12,:) == x(12,13,:) + x(12,14,:) + x(12,15,:) + x(12,17,:) + x(12,18,:) + x(12,19,:);
NETWORK.Constraints.flow10 = x(8,13,:) + x(12,13,:) >= x(13,16,:);
NETWORK.Constraints.flow11 = x(8,14,:) + x(12,14,:) >= x(14,16,:);
NETWORK.Constraints.flow12 = x(8,15,:) + x(12,15,:) >= x(15,16,:);
NETWORK.Constraints.flow13 = x(13,16,:) + x(14,16,:) + x(15,16,:) == x(16,17,:) + x(16,18,:) + x(16,19,:);
NETWORK.Constraints.flow14 = x(8,17,:) + x(12,17,:) + x(16,17,:) >= x(17,20,:);
NETWORK.Constraints.flow15 = x(8,18,:) + x(12,18,:) + x(16,18,:) >= x(18,20,:);
NETWORK.Constraints.flow16 = x(8,19,:) + x(12,19,:) + x(16,19,:) >= x(19,20,:);

%machine capacity
NETWORK.Constraints.machinecap1 = squeeze(x(4,5,:)) <= repmat(cm(1),[t,1]);
NETWORK.Constraints.machinecap2 = squeeze(x(4,6,:)) <= repmat(cm(2),[t,1]);
NETWORK.Constraints.machinecap3 = squeeze(x(4,7,:)) <= repmat(cm(3),[t,1]);
NETWORK.Constraints.machinecap4 = squeeze(x(8,9,:)) <= repmat(cm(4),[t,1]);
NETWORK.Constraints.machinecap5 = squeeze(x(8,10,:)) <= repmat(cm(5),[t,1]);
NETWORK.Constraints.machinecap6 = squeeze(x(8,11,:)) <= repmat(cm(6),[t,1]);
NETWORK.Constraints.machinecap7 = squeeze(x(8,13,:) + x(12,13,:)) <= repmat(cm(7),[t,1]);
NETWORK.Constraints.machinecap8 = squeeze(x(8,14,:) + x(12,14,:)) <= repmat(cm(8),[t,1]);
NETWORK.Constraints.machinecap9 = squeeze(x(8,15,:) + x(12,15,:)) <= repmat(cm(9),[t,1]);
NETWORK.Constraints.machinecap10 = squeeze(x(8,17,:) + x(12,17,:) + x(16,17,:)) <= repmat(cm(10),[t,1]);
NETWORK.Constraints.machinecap11 = squeeze(x(8,18,:) + x(12,18,:) + x(16,18,:)) <= repmat(cm(11),[t,1]);
NETWORK.Constraints.machinecap12 = squeeze(x(8,19,:) + x(12,19,:) + x(16,19,:)) <= repmat(cm(12),[t,1]);

% %storage capacity
NETWORK.Constraints.storagecapacity1 = squeeze(x(1,4,:) + x(2,4,:) + x(3,4,:)) <= repmat(cs(1),[t,1]);
NETWORK.Constraints.storagecapacity2 = squeeze(x(5,8,:) + x(6,8,:) + x(7,8,:)) <= repmat(cs(2),[t,1]);
NETWORK.Constraints.storagecapacity3 = squeeze(x(9,12,:) + x(10,12,:) + x(11,12,:)) <= repmat(cs(3),[t,1]);
NETWORK.Constraints.storagecapacity4 = squeeze(x(13,16,:) + x(14,16,:) + x(15,16,:)) <= repmat(cs(4),[t,1]);

NETWORK.ObjectiveSense = 'minimize';

NETWORK.Objective.COST = sum(x .* repmat(OCw,[1,1,t]),'all') + sum(x .* repmat(TCw,[1,1,t]),'all') + sum(a .* repmat(SC,[1,1,t]),'all') + sum(squeeze(sum(x,1)) .* repmat(HC,[1,t]),'all')...
                         + sum(R .* repmat(OCb,[1,1,t]),'all') + sum(S .* repmat(OCn,[1,1,t]),'all') + sum(V .* OCt,'all')...
                         + sum(R .* repmat(TCb,[1,1,t]),'all') + sum(S .* repmat(TCn,[1,1,t]),'all') + sum(T .* TCl,'all') + sum(U .* TCc,'all')...
                         + sum(R .* repmat(DCb,[1,1,t]),'all') + sum(T .* DCl,'all'); %+ sum(U .* SPc,'all')

[sol, fval, exitflag, output] = solve(NETWORK);


%%%%%%%%%%%%%%%%RATIO CONVERSIONS%%%%%%%%%%%%%%%%
flowratio = zeros(20,20,3);
totprelim = sol.x(4,5,:) + sol.x(4,6,:) + sol.x(4,7,:);
flowratio(4,5,:) = sol.x(4,5,:) ./ totprelim;
flowratio(4,6,:) = sol.x(4,6,:) ./ totprelim;
flowratio(4,7,:) = sol.x(4,7,:) ./ totprelim;

totstor1toprim = sol.x(8,9,:) + sol.x(8,10,:) + sol.x(8,11,:);
if totstor1toprim == 0
    flowratio(8,9,:) = 0;
    flowratio(8,10,:) = 0;
    flowratio(8,11,:) = 0;
else
    flowratio(8,9,:) = sol.x(8,9,:) ./ totstor1toprim;
    flowratio(8,10,:) = sol.x(8,10,:) ./ totstor1toprim;
    flowratio(8,11,:) = sol.x(8,11,:) ./ totstor1toprim;
end

totstor1tosec = sol.x(8,13,:) + sol.x(8,14,:) + sol.x(8,15,:);
if totstor1tosec == 0
    flowratio(8,13,:) = 0;
    flowratio(8,14,:) = 0;
    flowratio(8,15,:) = 0;
else
    flowratio(8,13,:) = sol.x(8,13,:) ./ totstor1tosec;
    flowratio(8,14,:) = sol.x(8,14,:) ./ totstor1tosec;
    flowratio(8,15,:) = sol.x(8,15,:) ./ totstor1tosec;
end

totstor1toter = sol.x(8,17,:) + sol.x(8,18,:) + sol.x(8,19,:);
if totstor1toter == 0
    flowratio(8,17,:) = 0;
    flowratio(8,18,:) = 0;
    flowratio(8,19,:) = 0;
else
    flowratio(8,17,:) = sol.x(8,17,:) ./ totstor1toter;
    flowratio(8,18,:) = sol.x(8,18,:) ./ totstor1toter;
    flowratio(8,19,:) = sol.x(8,19,:) ./ totstor1toter;
end

totstor2tosec = sol.x(12,13,:) + sol.x(12,14,:) + sol.x(12,15,:);
if totstor2tosec == 0
    flowratio(12,13,:) = 0;
    flowratio(12,14,:) = 0;
    flowratio(12,15,:) = 0;
else
    flowratio(12,13,:) = sol.x(12,13,:) ./ totstor2tosec;
    flowratio(12,14,:) = sol.x(12,14,:) ./ totstor2tosec;
    flowratio(12,15,:) = sol.x(12,15,:) ./ totstor2tosec;
end

totstor2toter = sol.x(12,17,:) + sol.x(12,18,:) + sol.x(12,19,:);
if totstor2toter == 0
    flowratio(12,17,:) = 0;
    flowratio(12,18,:) = 0;
    flowratio(12,19,:) = 0;
else
    flowratio(12,17,:) = sol.x(12,17,:) ./ totstor2toter;
    flowratio(12,18,:) = sol.x(12,18,:) ./ totstor2toter;
    flowratio(12,19,:) = sol.x(12,19,:) ./ totstor2toter;
end

totstor3toter = sol.x(16,17,:) + sol.x(16,18,:) + sol.x(16,19,:);
if totstor3toter == 0
    flowratio(16,17,:) = 0;
    flowratio(16,18,:) = 0;
    flowratio(16,19,:) = 0;
else
    flowratio(16,17,:) = sol.x(16,17,:) ./ totstor3toter;
    flowratio(16,18,:) = sol.x(16,18,:) ./ totstor3toter;
    flowratio(16,19,:) = sol.x(16,19,:) ./ totstor3toter;
end

WWTPCALC = optimproblem;

WWTPCALC.Constraints.trial1 = R .* repmat(domR,[1,1,t]) == RR;

%Byproduct disposal system
WWTPCALC.Constraints.BDS1 =    squeeze(sum(J,1)) == squeeze(sum(R,2));
WWTPCALC.Constraints.BDS2 =    squeeze(sum(R(:,1,:),1)) .* repmat(Q(1),[1,t])' == B(1,:)';
WWTPCALC.Constraints.BDS3 =    squeeze(sum(R(:,1,:),1)) .* repmat(Q(2),[1,t])' == B(2,:)';
WWTPCALC.Constraints.BDS4 =     B >= squeeze(sum(S,2));
WWTPCALC.Constraints.BDS5 =     squeeze(sum(R(:,2,:),1)) + squeeze(sum(S(:,2,:),1)) == sum(T,1)';
WWTPCALC.Constraints.BDS6 =     squeeze(sum(R(:,3,:),1)) + squeeze(sum(S(:,3,:),1)) == sum(U,1)';
WWTPCALC.Constraints.BDS7 =     T + udl >= dl;
WWTPCALC.Constraints.BDS8 =     U + udc >= dc;
WWTPCALC.Constraints.BDS9 =     squeeze(sum(R,1)) + squeeze(sum(S,1)) <= repmat(cp,[1,t]);
WWTPCALC.Constraints.BDS10 =     squeeze(sum(R(:,5,:),1)) + squeeze(sum(S(:,5,:),1)) == sum(V,1)';
WWTPCALC.Constraints.BDS11 =     A + udu >= du;
WWTPCALC.Constraints.BDS12 =     (squeeze(sum(S(:,4,:),1)) .* repmat(C,[1,t])') + (sum(V .* repmat(D,[1,t]) ,1)') == A;

%Effectiveness Calculator
%-------INITIAL COMPUTATION--------%
% computation of average input percentage for qp (aipq)
aipq(1,:) = sum(f .* squeeze(L(:,1,:)),1) ./ sum(f,1);
aipq(2,:) = sum(f .* squeeze(L(:,2,:)),1) ./ sum(f,1);
aipq(3,:) = sum(f .* squeeze(L(:,3,:)),1) ./ sum(f,1);
aipq(4,:) = sum(f .* squeeze(L(:,4,:)),1) ./ sum(f,1);

%computation of average input percentage for bp (aipb)
aipb(1,:) = sum(f .* squeeze(K(:,1,:)),1) ./ sum(f,1);
aipb(2,:) = sum(f .* squeeze(K(:,2,:)),1) ./ sum(f,1);
aipb(3,:) = sum(f .* squeeze(K(:,3,:)),1) ./ sum(f,1);

%-------INPUT--------%
WWTPCALC.Constraints.ina = f(1,:)' == squeeze(x(1,4,:));
WWTPCALC.Constraints.inb = f(2,:)' == squeeze(x(2,4,:));
WWTPCALC.Constraints.inc = f(3,:)' == squeeze(x(3,4,:));
WWTPCALC.Constraints.ind = (x(1,4,:) + x(2,4,:) + x(3,4,:)) .* flowratio(4,5,:) == x(4,5,:);
WWTPCALC.Constraints.ine = (x(1,4,:) + x(2,4,:) + x(3,4,:)) .* flowratio(4,6,:) == x(4,6,:);
WWTPCALC.Constraints.inf = (x(1,4,:) + x(2,4,:) + x(3,4,:)) .* flowratio(4,7,:) == x(4,7,:);

%-------FROM STORAGE 1--------%
prelimqpvol = optimconstr(j,l,t);
for j = 5:7
    for l = 1:l
        prelimqpvol(j,l,:) = squeeze(x(4,j,:)) .* aipq(l,:)' .* (1 - squeeze(repmat(Eq(j,l),[1,1,t]))) == squeeze(F(j,l,:));
    end
end
WWTPCALC.Constraints.fromstorage1 = prelimqpvol;

prelimbpvol = optimconstr(j,m,t);
for j = 5:7
    for m = 1:m
        prelimbpvol(j,m,:) = squeeze(x(4,j,:)) .* aipb(m,:)' .* (1 - squeeze(repmat(Eb(j,m),[1,1,t]))) == squeeze(G(j,m,:));
    end
end
WWTPCALC.Constraints.fromstorage2 = prelimbpvol;

prelimtobds = optimconstr(j,m,t);
for j = 5:7
    for m = 1:m
        prelimtobds(j,m,:) = squeeze(x(4,j,:)) .* aipb(m,:)' .* (squeeze(repmat(Eb(j,m),[1,1,t]))) == squeeze(J(j,m,:));
    end
end
WWTPCALC.Constraints.fromstorage3 = prelimtobds;

prelimclean = optimconstr(j,t);
for j = 5:7
    prelimclean(j,:) = squeeze(x(4,j,:)) .* (1 - sum(aipq,1)' - sum(aipb,1)') == H(j,:)';
end
WWTPCALC.Constraints.fromstorage4 = prelimclean;

stor2total = optimconstr(j,t);
for j = 5:7
    stor2total(j,:) = squeeze(x(j,8,:)) == squeeze(sum(F(j,:,:),2)) + squeeze(sum(G(j,:,:),2)) + H(j,:)';
end
WWTPCALC.Constraints.fromstorage5 = stor2total;

%-------FROM STORAGE 2--------%
frstor2dec = optimconstr(j,j,t);
for j = 9:11
    frstor2dec(8,j,:) = (x(5,8,:) + x(6,8,:) + x(7,8,:)) .* flowratio(8,j,:) == x(8,j,:);
end
for j = 13:15
    frstor2dec(8,j,:) = (x(5,8,:) + x(6,8,:) + x(7,8,:)) .* flowratio(8,j,:) == x(8,j,:);
end
for j = 17:19
    frstor2dec(8,j,:) = (x(5,8,:) + x(6,8,:) + x(7,8,:)) .* flowratio(8,j,:) == x(8,j,:);
end
WWTPCALC.Constraints.fromstorage6 = frstor2dec;

frstor2qpvol = optimconstr(j,l,t);
for l = 1:l
    for j = 9:11
        frstor2qpvol(j,l,:) = squeeze(F(5,l,:) + F(6,l,:) + F(7,l,:)) .* squeeze(flowratio(8,j,:)) .* (1 - squeeze(repmat(Eq(j,l),[1,1,t]))) == squeeze(F(j,l,:));
    end
    for j = 13:15
        frstor2qpvol(j,l,:) = squeeze(F(5,l,:) + F(6,l,:) + F(7,l,:)) .* squeeze(flowratio(8,j,:)) .* (1 - squeeze(repmat(Eq(j,l),[1,1,t]))) + squeeze(F(9,l,:) + F(10,l,:) + F(11,l,:)) .* squeeze(flowratio(12,j,:)) .* (1 - squeeze(repmat(Eq(j,l),[1,1,t]))) == squeeze(F(j,l,:));
    end
    for j = 17:19
        frstor2qpvol(j,l,:) = squeeze(F(5,l,:) + F(6,l,:) + F(7,l,:)) .* squeeze(flowratio(8,j,:)) .* (1 - squeeze(repmat(Eq(j,l),[1,1,t]))) + squeeze(F(9,l,:) + F(10,l,:) + F(11,l,:)) .* squeeze(flowratio(12,j,:)) .* (1 - squeeze(repmat(Eq(j,l),[1,1,t]))) + squeeze(F(13,l,:) + F(14,l,:) + F(15,l,:)) .* squeeze(flowratio(16,j,:)) .* (1 - squeeze(repmat(Eq(j,l),[1,1,t]))) == squeeze(F(j,l,:));
    end
end
WWTPCALC.Constraints.fromstorage7 = frstor2qpvol;

frstor2bpvol = optimconstr(j,m,t);
for m = 1:m
    for j = 9:11
        frstor2bpvol(j,m,:) = squeeze(F(5,m,:) + F(6,m,:) + F(7,m,:)) .* squeeze(flowratio(8,j,:)) .* (1 - squeeze(repmat(Eb(j,m),[1,1,t]))) == squeeze(G(j,m,:));
    end
    for j = 13:15
        frstor2bpvol(j,m,:) = squeeze(F(5,m,:) + F(6,m,:) + F(7,m,:)) .* squeeze(flowratio(8,j,:)) .* (1 - squeeze(repmat(Eb(j,m),[1,1,t]))) + squeeze(F(9,m,:) + F(10,m,:) + F(11,m,:)) .* squeeze(flowratio(12,j,:)) .* (1 - squeeze(repmat(Eb(j,m),[1,1,t]))) == squeeze(G(j,m,:));
    end
    for j = 17:19
        frstor2bpvol(j,m,:) = squeeze(F(5,m,:) + F(6,m,:) + F(7,m,:)) .* squeeze(flowratio(8,j,:)) .* (1 - squeeze(repmat(Eb(j,m),[1,1,t]))) + squeeze(F(9,m,:) + F(10,m,:) + F(11,m,:)) .* squeeze(flowratio(12,j,:)) .* (1 - squeeze(repmat(Eb(j,m),[1,1,t]))) + squeeze(F(13,m,:) + F(14,m,:) + F(15,m,:)) .* squeeze(flowratio(16,j,:)) .* (1 - squeeze(repmat(Eb(j,m),[1,1,t]))) == squeeze(G(j,m,:));
    end
end
WWTPCALC.Constraints.fromstorage8 = frstor2bpvol;

frstor2tobds = optimconstr(j,m,t);
for m = 1:m
    for j = 9:11
        frstor2tobds(j,m,:) = squeeze(F(5,m,:) + F(6,m,:) + F(7,m,:)) .* squeeze(flowratio(8,j,:)) .* (squeeze(repmat(Eb(j,m),[1,1,t]))) == squeeze(J(j,m,:));
    end
    for j = 13:15
        frstor2tobds(j,m,:) = squeeze(F(5,m,:) + F(6,m,:) + F(7,m,:)) .* squeeze(flowratio(8,j,:)) .* (squeeze(repmat(Eb(j,m),[1,1,t]))) + squeeze(F(9,m,:) + F(10,m,:) + F(11,m,:)) .* squeeze(flowratio(12,j,:)) .* (squeeze(repmat(Eb(j,m),[1,1,t]))) == squeeze(J(j,m,:));
    end
    for j = 17:19
        frstor2tobds(j,m,:) = squeeze(F(5,m,:) + F(6,m,:) + F(7,m,:)) .* squeeze(flowratio(8,j,:)) .* (squeeze(repmat(Eb(j,m),[1,1,t]))) + squeeze(F(9,m,:) + F(10,m,:) + F(11,m,:)) .* squeeze(flowratio(12,j,:)) .* (squeeze(repmat(Eb(j,m),[1,1,t]))) + squeeze(F(13,m,:) + F(14,m,:) + F(15,m,:)) .* squeeze(flowratio(16,j,:)) .* (squeeze(repmat(Eb(j,m),[1,1,t])))== squeeze(J(j,m,:));
    end
end
WWTPCALC.Constraints.fromstorage9 = frstor2tobds;

primclean = optimconstr(j,t);
for j = 9:11
    primclean(j,:) = squeeze(x(8,j,:)) - squeeze(flowratio(8,j,:)) .* squeeze((sum(F(5,:,:),2) + sum(F(6,:,:),2) + sum(F(7,:,:),2) + sum(G(5,:,:),2) + sum(G(6,:,:),2) + sum(G(7,:,:),2))) == H(j,:)';
end
WWTPCALC.Constraints.fromstorage10 = primclean;

stor3total = optimconstr(j,t);
for j = 9:11
    stor3total(j,:) = squeeze(x(j,12,:)) == squeeze(sum(F(j,:,:),2)) + squeeze(sum(G(j,:,:),2)) + H(j,:)';
end
WWTPCALC.Constraints.fromstorage11 = stor3total;


%-------FROM STORAGE 3--------%
frstor3dec = optimconstr(j,j,t);
for j = 13:15
    frstor3dec(12,j,:) = (x(9,12,:) + x(10,12,:) + x(11,12,:)) .* flowratio(12,j,:) == x(12,j,:);
end
for j = 17:19
    frstor3dec(12,j,:) = (x(9,12,:) + x(10,12,:) + x(11,12,:)) .* flowratio(12,j,:) == x(12,j,:);
end
WWTPCALC.Constraints.fromstorage12 = frstor3dec;

secclean = optimconstr(j,t);
for j = 13:15
    secclean(j,:) = (squeeze(x(8,j,:)) - squeeze(flowratio(8,j,:)) .* squeeze((sum(F(5,:,:),2) + sum(F(6,:,:),2) + sum(F(7,:,:),2) + sum(G(5,:,:),2) + sum(G(6,:,:),2) + sum(G(7,:,:),2)))) + (squeeze(x(12,j,:)) - squeeze(flowratio(12,j,:)) .* squeeze((sum(F(9,:,:),2) + sum(F(10,:,:),2) + sum(F(11,:,:),2) + sum(G(9,:,:),2) + sum(G(10,:,:),2) + sum(G(11,:,:),2)))) == H(j,:)';
end
WWTPCALC.Constraints.fromstorage13 = secclean;

stor4total = optimconstr(j,t);
for j = 13:15
    stor4total(j,:) = squeeze(x(j,16,:)) == squeeze(sum(F(j,:,:),2)) + squeeze(sum(G(j,:,:),2)) + H(j,:)';
end
WWTPCALC.Constraints.fromstorage14 = stor4total;

%-------FROM STORAGE 4--------%
frstor4dec = optimconstr(j,j,t);
for j = 17:19
    frstor4dec(16,j,:) = (x(13,16,:) + x(14,16,:) + x(15,16,:)) .* flowratio(16,j,:) == x(16,j,:);
end
WWTPCALC.Constraints.fromstorage15 = frstor4dec;

terclean = optimconstr(j,t);
for j = 17:19
    terclean(j,:) = (squeeze(x(8,j,:)) - squeeze(flowratio(8,j,:)) .* squeeze((sum(F(5,:,:),2) + sum(F(6,:,:),2) + sum(F(7,:,:),2) + sum(G(5,:,:),2) + sum(G(6,:,:),2) + sum(G(7,:,:),2)))) + (squeeze(x(12,j,:)) - squeeze(flowratio(12,j,:)) .* squeeze((sum(F(9,:,:),2) + sum(F(10,:,:),2) + sum(F(11,:,:),2) + sum(G(9,:,:),2) + sum(G(10,:,:),2) + sum(G(11,:,:),2)))) + (squeeze(x(16,j,:)) - squeeze(flowratio(16,j,:)) .* squeeze((sum(F(13,:,:),2) + sum(F(14,:,:),2) + sum(F(15,:,:),2) + sum(G(13,:,:),2) + sum(G(14,:,:),2) + sum(G(15,:,:),2)))) == H(j,:)';
end
WWTPCALC.Constraints.fromstorage16 = terclean;

outputtotal = optimconstr(j,t);
for j = 17:19
    outputtotal(j,:) = squeeze(x(j,20,:)) == squeeze(sum(F(j,:,:),2)) + squeeze(sum(G(j,:,:),2)) + H(j,:)';
end
WWTPCALC.Constraints.fromstorage17 = outputtotal;

WWTPCALC.Constraints.MULTOBJVAR = ENVI == sum(x .* repmat(PEw,[1,1,t]),'all') + sum(x .* repmat(TEw,[1,1,t]),'all')...
                                    + sum(R .* repmat(PEb,[1,1,t]),'all') + sum(R .* repmat(DEb,[1,1,t]),'all')...
                                    + sum(S .* repmat(PEn,[1,1,t]),'all') + sum(S .* repmat(DEn,[1,1,t]),'all')...
                                    + sum(V .* repmat(PEt,[1,t]),'all') + sum(U .* repmat(PEc,[1,t]),'all') + sum(T .* repmat(DEl,[1,t]),'all')...
                                    - sum(S .* repmat(GHG,[1,1,t]),'all');

WWTPCALC.ObjectiveSense = 'minimize';

WWTPCALC.Objective.COST = sum(x .* repmat(OCw,[1,1,t]),'all') + sum(x .* repmat(TCw,[1,1,t]),'all') + sum(squeeze(sum(x,1)) .* repmat(HC,[1,t]),'all')...
                         + sum(R .* repmat(OCb,[1,1,t]),'all') + sum(S .* repmat(OCn,[1,1,t]),'all') + sum(V .* OCt,'all')...
                         + sum(R .* repmat(TCb,[1,1,t]),'all') + sum(S .* repmat(TCn,[1,1,t]),'all') + sum(T .* TCl,'all') + sum(U .* TCc,'all')...
                         + sum(R .* repmat(DCb,[1,1,t]),'all') + sum(T .* DCl,'all') - sum(U .* SPc,'all')...
                         + ENVI*0.3;

[sol1, fval1, exitflag1, output1] = solve(WWTPCALC);
fval1

ENVPEw = sum(sol1.x .* repmat(PEw,[1,1,t]),'all');
ENVTEw = sum(sol1.x .* repmat(TEw,[1,1,t]),'all');
ENVPEr = sum(sol1.R .* repmat(PEb,[1,1,t]),'all');
ENVPEs = sum(sol1.S .* repmat(PEn,[1,1,t]),'all');
ENVPEv = sum(sol1.V .* repmat(PEt,[1,t]),'all');
ENVDEr =  sum(sol1.R .* repmat(DEb,[1,1,t]),'all');
ENVDEs = sum(sol1.S .* repmat(DEn,[1,1,t]),'all');
ENVPEu = sum(sol1.U .* repmat(PEc,[1,t]),'all');
ENVDEt = sum(sol1.T .* repmat(DEl,[1,t]),'all');
ENVSave = sum(sol1.S .* repmat(GHG,[1,1,t]),'all');

NCALC = zeros(20,4,3);
for l = 1:l
    %node 8
    if sol1.x(5,8,:) + sol1.x(6,8,:) + sol1.x(7,8,:) == 0
        NCALC(8,l,:) = 0;
    else
        NCALC(8,l,:) = (sol1.F(5,l,:) + sol1.F(6,l,:) + sol1.F(7,l,:) ) ./ (sol1.x(5,8,:) + sol1.x(6,8,:) + sol1.x(7,8,:));
    end
    %node 12
    if sol1.x(9,12,:) + sol1.x(10,12,:) + sol1.x(11,12,:) == 0
        NCALC(12,l,:) = 0;
    else
        NCALC(12,l,:) = (sol1.F(9,l,:) + sol1.F(10,l,:) + sol1.F(11,l,:) ) ./ (sol1.x(9,12,:) + sol1.x(10,12,:) + sol1.x(11,12,:));
    end
    %node 16
    if sol1.x(13,16,:) + sol1.x(14,16,:) + sol1.x(15,16,:) == 0
        NCALC(16,l,:) = 0;
    else
        NCALC(16,l,:) = (sol1.F(13,l,:) + sol1.F(14,l,:) + sol1.F(15,l,:) ) ./ (sol1.x(13,16,:) + sol1.x(14,16,:) + sol1.x(15,16,:));
    end
    %node 20
    if sol1.x(17,20,:) + sol1.x(18,20,:) + sol1.x(19,20,:) == 0
        NCALC(20,l,:) = 0;
    else
        NCALC(20,l,:) = (sol1.F(17,l,:) + sol1.F(18,l,:) + sol1.F(19,l,:) ) ./ (sol1.x(17,20,:) + sol1.x(18,20,:) + sol1.x(19,20,:));
    end
end

OCALC = zeros(20,3,3);
for m = 1:m
    %node 8
    if sol1.x(5,8,:) + sol1.x(6,8,:) + sol1.x(7,8,:) == 0
        OCALC(8,m,:) = 0;
    else
        OCALC(8,m,:) = (sol1.G(5,m,:) + sol1.G(6,m,:) + sol1.G(7,m,:) ) ./ (sol1.x(5,8,:) + sol1.x(6,8,:) + sol1.x(7,8,:));
    end
    %node 12
    if sol1.x(9,12,:) + sol1.x(10,12,:) + sol1.x(11,12,:) == 0
        OCALC(12,m,:) = 0;
    else
        OCALC(12,m,:) = (sol1.G(9,m,:) + sol1.G(10,m,:) + sol1.G(11,m,:) ) ./ (sol1.x(9,12,:) + sol1.x(10,12,:) + sol1.x(11,12,:));
    end
    %node 16
    if sol1.x(13,16,:) + sol1.x(14,16,:) + sol1.x(15,16,:) == 0
        OCALC(16,m,:) = 0;
    else
        OCALC(16,m,:) = (sol1.G(13,m,:) + sol1.G(14,m,:) + sol1.G(15,m,:) ) ./ (sol1.x(13,16,:) + sol1.x(14,16,:) + sol1.x(15,16,:));
    end
    %node 20
    if sol1.x(17,20,:) + sol1.x(18,20,:) + sol1.x(19,20,:) == 0
        OCALC(20,m,:) = 0;
    else
        OCALC(20,m,:) = (sol1.G(17,m,:) + sol1.G(18,m,:) + sol1.G(19,m,:) ) ./ (sol1.x(17,20,:) + sol1.x(18,20,:) + sol1.x(19,20,:));
    end
end


i = 1;
[sol2, fval2, exitflag2, output2] = N_O_parameter_function(NCALC,OCALC,Eq,Eb);

xdum = zeros(20,20,3);
xdum(:,:,:,i) = sol2.x(:,:,:,i);
xdum1 = squeeze(xdum(:,:,:,i));
xdum2 = zeros(20,20,3);
Ndum = zeros(20,4,3);
Odum = zeros(20,3,3);

while true

    [sol3, fval3, exitflag3, output3] = x_parameter_function(xdum1,Eq,Eb);
    Ndum(:,:,:,i) = sol3.N;
    Ndum1 = squeeze(Ndum(:,:,:,i));
    Odum(:,:,:,i) = sol3.O;
    Odum1 = squeeze(Odum(:,:,:,i));
    [sol4, fval4, exitflag4, output4] = N_O_parameter_function(Ndum1,Odum1,Eq,Eb);
    xdum2(:,:,:,i) = sol4.x;
    xdum3 = squeeze(xdum2(:,:,:,i));

    if abs(xdum3 - xdum1) <= 0.5
        break;
    else
        xdum1 = xdum3;
    end

    i=i+1;
end

fx = sol4.x;
fa = sol4.a;
fb = sol4.b;
fe = sol4.e;
fR = sol4.R;
fS = sol4.S;
fT = sol4.T;
fU = sol4.U;
fV = sol4.V;
fF = sol4.F;
fG = sol4.G;
fH = sol4.H;
fJ = sol4.J;
fN = sol3.N;
fO = sol3.O;
fB = sol4.B;
fA = sol4.A;
fENVI = sol4.ENVI;

optimsol = sum(fx .* repmat(OCw,[1,1,t]),'all') + sum(fx .* repmat(TCw,[1,1,t]),'all') + sum(fa .* repmat(SC,[1,1,t]),'all') + sum(squeeze(sum(fx,2)) .* repmat(HC,[1,t]),'all')... 
         + sum(fR .* repmat(OCb,[1,1,t]),'all') + sum(fS .* repmat(OCn,[1,1,t]),'all') + sum(fV .* OCt,'all')...
         + sum(fR .* repmat(TCb,[1,1,t]),'all') + sum(fS .* repmat(TCn,[1,1,t]),'all') + sum(fT .* TCl,'all') + sum(fU .* TCc,'all')...
         + sum(fR .* repmat(DCb,[1,1,t]),'all') + sum(fT .* DCl,'all') - sum(fU .* SPc,'all')...
         + fENVI*0.3;

TOCW = sum(fx .* repmat(OCw,[1,1,t]),'all');
TOCb = sum(fR .* repmat(OCb,[1,1,t]),'all') + sum(fS .* repmat(OCn,[1,1,t]),'all') + sum(fV .* OCt,'all');
TTCw = sum(fx .* repmat(TCw,[1,1,t]),'all');
TSCq = sum(fa .* repmat(SC,[1,1,t]),'all');
THCw = sum(squeeze(sum(fx,2)) .* repmat(HC,[1,t]),'all');
TTCb = sum(fR .* repmat(TCb,[1,1,t]),'all') + sum(fS .* repmat(TCn,[1,1,t]),'all') + sum(fT .* TCl,'all') + sum(fU .* TCc,'all');
TDCb = sum(fR .* repmat(DCb,[1,1,t]),'all') + sum(fT .* DCl,'all');
SPu = sum(fU .* SPc,'all');
PENALTY = fENVI*0.3;

writematrix(fx(:,:,1),'x.xlsx','Sheet','Sheet1','Range','A1:T20');
writematrix(fx(:,:,2),'x.xlsx','Sheet','Sheet1','Range','A23:T42');
writematrix(fx(:,:,3),'x.xlsx','Sheet','Sheet1','Range','A45:T64'); 


%results analyzer
writematrix(optimsol,'results_analyzer.xlsx','Sheet','ObjFunc','Range','B1');
writematrix(fENVI,'results_analyzer.xlsx','Sheet','ObjFunc','Range','B2');
writematrix("Optimal Solution",'results_analyzer.xlsx','Sheet','ObjFunc','Range','A1');
writematrix("Environmental Emissions",'results_analyzer.xlsx','Sheet','ObjFunc','Range','A2');

%sol.x
writematrix(["Period 1" "Machine 1" "Machine 2" "Machine 3"], 'results_analyzer.xlsx', 'Sheet','sol.x','Range','A2:D2');
writematrix(["Period 2" "Machine 1" "Machine 2" "Machine 3"], 'results_analyzer.xlsx', 'Sheet','sol.x','Range','A14:D14');
writematrix(["Period 3" "Machine 1" "Machine 2" "Machine 3"], 'results_analyzer.xlsx', 'Sheet','sol.x','Range','A26:D26');
writematrix(["Storage 1"; "Preliminary"; "Storage 2"; "Primary"; "Storage 3"; "Secondary"; "Storage 4"; "Tertiary"], 'results_analyzer.xlsx', 'Sheet','sol.x','Range','A3:A10');
writematrix(["Storage 1"; "Preliminary"; "Storage 2"; "Primary"; "Storage 3"; "Secondary"; "Storage 4"; "Tertiary"], 'results_analyzer.xlsx', 'Sheet','sol.x','Range','A15:A22');
writematrix(["Storage 1"; "Preliminary"; "Storage 2"; "Primary"; "Storage 3"; "Secondary"; "Storage 4"; "Tertiary"], 'results_analyzer.xlsx', 'Sheet','sol.x','Range','A27:A34');

writematrix([fx(4,5,1) fx(4,6,1) fx(4,7,1)],'results_analyzer.xlsx', 'Sheet','sol.x','Range','B3:D3');
writematrix([fx(5,8,1) fx(6,8,1) fx(7,8,1)],'results_analyzer.xlsx', 'Sheet','sol.x','Range','B4:D4');
writematrix([fx(8,9,1) fx(8,10,1) fx(8,11,1)],'results_analyzer.xlsx', 'Sheet','sol.x','Range','B5:D5');
writematrix([fx(9,12,1) fx(10,12,1) fx(11,12,1)],'results_analyzer.xlsx', 'Sheet','sol.x','Range','B6:D6');
writematrix([fx(12,13,1) fx(12,14,1) fx(12,15,1)],'results_analyzer.xlsx', 'Sheet','sol.x','Range','B7:D7');
writematrix([fx(13,16,1) fx(14,16,1) fx(15,16,1)],'results_analyzer.xlsx', 'Sheet','sol.x','Range','B8:D8');
writematrix([fx(16,17,1) fx(16,18,1) fx(16,19,1)],'results_analyzer.xlsx', 'Sheet','sol.x','Range','B9:D9');
writematrix([fx(17,20,1) fx(18,20,1) fx(19,20,1)],'results_analyzer.xlsx', 'Sheet','sol.x','Range','B10:D10');

writematrix([fx(4,5,2) fx(4,6,2) fx(4,7,2)],'results_analyzer.xlsx', 'Sheet','sol.x','Range','B15:D15');
writematrix([fx(5,8,2) fx(6,8,2) fx(7,8,2)],'results_analyzer.xlsx', 'Sheet','sol.x','Range','B16:D16');
writematrix([fx(8,9,2) fx(8,10,2) fx(8,11,2)],'results_analyzer.xlsx', 'Sheet','sol.x','Range','B17:D17');
writematrix([fx(9,12,2) fx(10,12,2) fx(11,12,2)],'results_analyzer.xlsx', 'Sheet','sol.x','Range','B18:D18');
writematrix([fx(12,13,2) fx(12,14,2) fx(12,15,2)],'results_analyzer.xlsx', 'Sheet','sol.x','Range','B19:D19');
writematrix([fx(13,16,2) fx(14,16,2) fx(15,16,2)],'results_analyzer.xlsx', 'Sheet','sol.x','Range','B20:D20');
writematrix([fx(16,17,2) fx(16,18,2) fx(16,19,2)],'results_analyzer.xlsx', 'Sheet','sol.x','Range','B21:D21');
writematrix([fx(17,20,2) fx(18,20,2) fx(19,20,2)],'results_analyzer.xlsx', 'Sheet','sol.x','Range','B22:D22');

writematrix([fx(4,5,3) fx(4,6,3) fx(4,7,3)],'results_analyzer.xlsx', 'Sheet','sol.x','Range','B27:D27');
writematrix([fx(5,8,3) fx(6,8,3) fx(7,8,3)],'results_analyzer.xlsx', 'Sheet','sol.x','Range','B28:D28');
writematrix([fx(8,9,3) fx(8,10,3) fx(8,11,3)],'results_analyzer.xlsx', 'Sheet','sol.x','Range','B29:D29');
writematrix([fx(9,12,3) fx(10,12,3) fx(11,12,3)],'results_analyzer.xlsx', 'Sheet','sol.x','Range','B30:D30');
writematrix([fx(12,13,3) fx(12,14,3) fx(12,15,3)],'results_analyzer.xlsx', 'Sheet','sol.x','Range','B31:D31');
writematrix([fx(13,16,3) fx(14,16,3) fx(15,16,3)],'results_analyzer.xlsx', 'Sheet','sol.x','Range','B32:D32');
writematrix([fx(16,17,3) fx(16,18,3) fx(16,19,3)],'results_analyzer.xlsx', 'Sheet','sol.x','Range','B33:D33');
writematrix([fx(17,20,3) fx(18,20,3) fx(19,20,3)],'results_analyzer.xlsx', 'Sheet','sol.x','Range','B34:D34');

%sol.R
writematrix(["Period 1" "ADP" "L" "C" "BU" "I"], 'results_analyzer.xlsx', 'Sheet','sol.R','Range','A2:F2');
writematrix(["LSC"; "S"; "AS";], 'results_analyzer.xlsx', 'Sheet','sol.R','Range','A3:A5');
writematrix([fR(1,1,1) fR(1,2,1) fR(1,3,1) fR(1,4,1) fR(1,5,1)],'results_analyzer.xlsx', 'Sheet','sol.R','Range','B3:F3');
writematrix([fR(2,1,1) fR(2,2,1) fR(2,3,1) fR(2,4,1) fR(2,5,1)],'results_analyzer.xlsx', 'Sheet','sol.R','Range','B4:F4');
writematrix([fR(3,1,1) fR(3,2,1) fR(3,3,1) fR(3,4,1) fR(3,5,1)],'results_analyzer.xlsx', 'Sheet','sol.R','Range','B5:F5');

writematrix(["Period 2" "ADP" "L" "C" "BU" "I"], 'results_analyzer.xlsx', 'Sheet','sol.R','Range','A10:F10');
writematrix(["LSC"; "S"; "AS";], 'results_analyzer.xlsx', 'Sheet','sol.R','Range','A11:A13');
writematrix([fR(1,1,2) fR(1,2,2) fR(1,3,2) fR(1,4,2) fR(1,5,2)],'results_analyzer.xlsx', 'Sheet','sol.R','Range','B11:F11');
writematrix([fR(2,1,2) fR(2,2,2) fR(2,3,2) fR(2,4,2) fR(2,5,2)],'results_analyzer.xlsx', 'Sheet','sol.R','Range','B12:F12');
writematrix([fR(3,1,2) fR(3,2,2) fR(3,3,2) fR(3,4,2) fR(3,5,2)],'results_analyzer.xlsx', 'Sheet','sol.R','Range','B13:F13');


writematrix(["Period 3" "ADP" "L" "C" "BU" "I"], 'results_analyzer.xlsx', 'Sheet','sol.R','Range','A17:F17');
writematrix(["LSC"; "S"; "AS";], 'results_analyzer.xlsx', 'Sheet','sol.R','Range','A18:A20');
writematrix([fR(1,1,3) fR(1,2,3) fR(1,3,3) fR(1,4,3) fR(1,5,3)],'results_analyzer.xlsx', 'Sheet','sol.R','Range','B18:F18');
writematrix([fR(2,1,3) fR(2,2,3) fR(2,3,3) fR(2,4,3) fR(2,5,3)],'results_analyzer.xlsx', 'Sheet','sol.R','Range','B19:F19');
writematrix([fR(3,1,3) fR(3,2,3) fR(3,3,3) fR(3,4,3) fR(3,5,3)],'results_analyzer.xlsx', 'Sheet','sol.R','Range','B20:F20');

%sol.S
writematrix(["Period 1" "ADP" "L" "C" "BU" "I"], 'results_analyzer.xlsx', 'Sheet','sol.S','Range','A2:F2');
writematrix(["Biogas"; "Digestate";], 'results_analyzer.xlsx', 'Sheet','sol.S','Range','A3:A4');
writematrix([fS(1,1,1) fS(1,2,1) fS(1,3,1) fS(1,4,1) fS(1,5,1)],'results_analyzer.xlsx', 'Sheet','sol.S','Range','B3:F3');
writematrix([fS(2,1,1) fS(2,2,1) fS(2,3,1) fS(2,4,1) fS(2,5,1)],'results_analyzer.xlsx', 'Sheet','sol.S','Range','B4:F4');

writematrix(["Period 2" "ADP" "L" "C" "BU" "I"], 'results_analyzer.xlsx', 'Sheet','sol.S','Range','A8:F8');
writematrix(["Biogas"; "Digestate";], 'results_analyzer.xlsx', 'Sheet','sol.S','Range','A9:A10');
writematrix([fS(1,1,2) fS(1,2,2) fS(1,3,2) fS(1,4,2) fS(1,5,2)],'results_analyzer.xlsx', 'Sheet','sol.S','Range','B9:F9');
writematrix([fS(2,1,2) fS(2,2,2) fS(2,3,2) fS(2,4,2) fS(2,5,2)],'results_analyzer.xlsx', 'Sheet','sol.S','Range','B10:F10');


writematrix(["Period 3" "ADP" "L" "C" "BU" "I"], 'results_analyzer.xlsx', 'Sheet','sol.S','Range','A14:F14');
writematrix(["Biogas"; "Digestate";], 'results_analyzer.xlsx', 'Sheet','sol.S','Range','A15:A16');
writematrix([fS(1,1,3) fS(1,2,3) fS(1,3,3) fS(1,4,3) fS(1,5,3)],'results_analyzer.xlsx', 'Sheet','sol.S','Range','B15:F15');
writematrix([fS(2,1,3) fS(2,2,3) fS(2,3,3) fS(2,4,3) fS(2,5,3)],'results_analyzer.xlsx', 'Sheet','sol.S','Range','B16:F16');


%sol.T
writematrix(["Period 1" "Period 2" "Period 3"], 'results_analyzer.xlsx', 'Sheet','sol.T','Range','B1:F1');
writematrix(["Landfill 1"; "Landfill 2"; "Landfill 3"], 'results_analyzer.xlsx', 'Sheet','sol.T','Range','A2:A4');
writematrix([fT(1,1) fT(1,2) fT(1,3)],'results_analyzer.xlsx', 'Sheet','sol.T','Range','B2:D2');
writematrix([fT(2,1) fT(2,2) fT(2,3)],'results_analyzer.xlsx', 'Sheet','sol.T','Range','B3:D3');
writematrix([fT(3,1) fT(3,2) fT(3,3)],'results_analyzer.xlsx', 'Sheet','sol.T','Range','B4:D4');

%sol.U
writematrix(["Period 1" "Period 2" "Period 3"], 'results_analyzer.xlsx', 'Sheet','sol.U','Range','B1:F1');
writematrix(["Agricultural Buyer 1"; "Agricultural Buyer 2"; "Agricultural Buyer 3"], 'results_analyzer.xlsx', 'Sheet','sol.U','Range','A2:A4');
writematrix([fU(1,1) fU(1,2) fU(1,3)],'results_analyzer.xlsx', 'Sheet','sol.U','Range','B2:D2');
writematrix([fU(2,1) fU(2,2) fU(2,3)],'results_analyzer.xlsx', 'Sheet','sol.U','Range','B3:D3');
writematrix([fU(3,1) fU(3,2) fU(3,3)],'results_analyzer.xlsx', 'Sheet','sol.U','Range','B4:D4');

%sol.V
writematrix(["Period 1" "Period 2" "Period 3"], 'results_analyzer.xlsx', 'Sheet','sol.V','Range','B1:F1');
writematrix(["Turbine Generator"; "Thermoelectric Generator";], 'results_analyzer.xlsx', 'Sheet','sol.V','Range','A2:A3');
writematrix([fV(1,1) fV(1,2) fV(1,3)],'results_analyzer.xlsx', 'Sheet','sol.V','Range','B2:D2');
writematrix([fV(2,1) fV(2,2) fV(2,3)],'results_analyzer.xlsx', 'Sheet','sol.V','Range','B3:D3');

%sol.F
writematrix(["Period 1" "BOD" "COD" "N" "P"], 'results_analyzer.xlsx', 'Sheet','sol.F','Range','A1:E1');
writematrix(["Preliminary 1"; "Preliminary 2"; "Preliminary 3"; "Primary 1"; "Primary 2"; "Primary 3"; "Secondary 1"; "Secondary 2"; "Secondary 3"; "Tertiary 1"; "Tertiary 2"; "Tertiary 3" ], 'results_analyzer.xlsx', 'Sheet','sol.F','Range','A2:A13');
writematrix([fF(5,1,1) fF(5,2,1) fF(5,3,1) fF(5,4,1)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B2:E3');
writematrix([fF(6,1,1) fF(6,2,1) fF(6,3,1) fF(6,4,1)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B3:E3');
writematrix([fF(7,1,1) fF(7,2,1) fF(7,3,1) fF(7,4,1)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B4:E4');
writematrix([fF(9,1,1) fF(9,2,1) fF(9,3,1) fF(9,4,1)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B5:E5');
writematrix([fF(10,1,1) fF(10,2,1) fF(10,3,1) fF(10,4,1)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B6:E6');
writematrix([fF(11,1,1) fF(11,2,1) fF(11,3,1) fF(11,4,1)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B7:E7');
writematrix([fF(13,1,1) fF(13,2,1) fF(13,3,1) fF(13,4,1)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B8:E8');
writematrix([fF(14,1,1) fF(14,2,1) fF(14,3,1) fF(14,4,1)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B9:E9');
writematrix([fF(15,1,1) fF(15,2,1) fF(15,3,1) fF(15,4,1)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B10:E10');
writematrix([fF(17,1,1) fF(17,2,1) fF(17,3,1) fF(17,4,1)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B11:E11');
writematrix([fF(18,1,1) fF(18,2,1) fF(18,3,1) fF(18,4,1)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B12:E12');
writematrix([fF(19,1,1) fF(19,2,1) fF(19,3,1) fF(19,4,1)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B13:E13');


writematrix(["Period 2" "BOD" "COD" "N" "P"], 'results_analyzer.xlsx', 'Sheet','sol.F','Range','A17:E17');
writematrix(["Preliminary 1"; "Preliminary 2"; "Preliminary 3"; "Primary 1"; "Primary 2"; "Primary 3"; "Secondary 1"; "Secondary 2"; "Secondary 3"; "Tertiary 1"; "Tertiary 2"; "Tertiary 3" ], 'results_analyzer.xlsx', 'Sheet','sol.F','Range','A18:A29');
writematrix([fF(5,1,2) fF(5,2,2) fF(5,3,2) fF(5,4,2)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B18:E18');
writematrix([fF(6,1,2) fF(6,2,2) fF(6,3,2) fF(6,4,2)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B19:E19');
writematrix([fF(7,1,2) fF(7,2,2) fF(7,3,2) fF(7,4,2)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B20:E20');
writematrix([fF(9,1,2) fF(9,2,2) fF(9,3,2) fF(9,4,2)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B21:E21');
writematrix([fF(10,1,2) fF(10,2,2) fF(10,3,2) fF(10,4,2)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B22:E22');
writematrix([fF(11,1,2) fF(11,2,2) fF(11,3,2) fF(11,4,2)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B23:E23');
writematrix([fF(13,1,2) fF(13,2,2) fF(13,3,2) fF(13,4,2)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B24:E24');
writematrix([fF(14,1,2) fF(14,2,2) fF(14,3,2) fF(14,4,2)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B25:E25');
writematrix([fF(15,1,2) fF(15,2,2) fF(15,3,2) fF(15,4,2)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B26:E26');
writematrix([fF(17,1,2) fF(17,2,2) fF(17,3,2) fF(17,4,2)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B27:E27');
writematrix([fF(18,1,2) fF(18,2,2) fF(18,3,2) fF(18,4,2)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B28:E28');
writematrix([fF(19,1,2) fF(19,2,2) fF(19,3,2) fF(19,4,2)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B29:E29');


writematrix(["Period 3" "BOD" "COD" "N" "P"], 'results_analyzer.xlsx', 'Sheet','sol.F','Range','A33:E33');
writematrix(["Preliminary 1"; "Preliminary 2"; "Preliminary 3"; "Primary 1"; "Primary 2"; "Primary 3"; "Secondary 1"; "Secondary 2"; "Secondary 3"; "Tertiary 1"; "Tertiary 2"; "Tertiary 3" ], 'results_analyzer.xlsx', 'Sheet','sol.F','Range','A34:A45');
writematrix([fF(5,1,3) fF(5,2,3) fF(5,3,3) fF(5,4,3)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B34:E34');
writematrix([fF(6,1,3) fF(6,2,3) fF(6,3,3) fF(6,4,3)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B35:E35');
writematrix([fF(7,1,3) fF(7,2,3) fF(7,3,3) fF(7,4,3)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B36:E36');
writematrix([fF(9,1,3) fF(9,2,3) fF(9,3,3) fF(9,4,3)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B37:E37');
writematrix([fF(10,1,3) fF(10,2,3) fF(10,3,3) fF(10,4,3)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B38:E38');
writematrix([fF(11,1,3) fF(11,2,3) fF(11,3,3) fF(11,4,3)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B39:E39');
writematrix([fF(13,1,3) fF(13,2,3) fF(13,3,3) fF(13,4,3)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B40:E40');
writematrix([fF(14,1,3) fF(14,2,3) fF(14,3,3) fF(14,4,3)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B41:E41');
writematrix([fF(15,1,3) fF(15,2,3) fF(15,3,3) fF(15,4,3)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B42:E42');
writematrix([fF(17,1,3) fF(17,2,3) fF(17,3,3) fF(17,4,3)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B43:E43');
writematrix([fF(18,1,3) fF(18,2,3) fF(18,3,3) fF(18,4,3)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B44:E44');
writematrix([fF(19,1,3) fF(19,2,3) fF(19,3,3) fF(19,4,3)],'results_analyzer.xlsx', 'Sheet','sol.F','Range','B45:E45');

%sol.G
writematrix(["Period 1" "LSC" "Sludge" "Activated Sludge"], 'results_analyzer.xlsx', 'Sheet','sol.G','Range','A1:D1');
writematrix(["Preliminary 1"; "Preliminary 2"; "Preliminary 3"; "Primary 1"; "Primary 2"; "Primary 3"; "Secondary 1"; "Secondary 2"; "Secondary 3"; "Tertiary 1"; "Tertiary 2"; "Tertiary 3" ], 'results_analyzer.xlsx', 'Sheet','sol.G','Range','A2:A13');
writematrix([fG(5,1,1) fG(5,2,1) fG(5,3,1)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B2:D3');
writematrix([fG(6,1,1) fG(6,2,1) fG(6,3,1)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B3:D3');
writematrix([fG(7,1,1) fG(7,2,1) fG(7,3,1)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B4:D4');
writematrix([fG(9,1,1) fG(9,2,1) fG(9,3,1)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B5:D5');
writematrix([fG(10,1,1) fG(10,2,1) fG(10,3,1)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B6:D6');
writematrix([fG(11,1,1) fG(11,2,1) fG(11,3,1)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B7:D7');
writematrix([fG(13,1,1) fG(13,2,1) fG(13,3,1)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B8:D8');
writematrix([fG(14,1,1) fG(14,2,1) fG(14,3,1)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B9:D9');
writematrix([fG(15,1,1) fG(15,2,1) fG(15,3,1)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B10:D10');
writematrix([fG(17,1,1) fG(17,2,1) fG(17,3,1)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B11:D11');
writematrix([fG(18,1,1) fG(18,2,1) fG(18,3,1)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B12:D12');
writematrix([fG(19,1,1) fG(19,2,1) fG(19,3,1)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B13:D13');

writematrix(["Period 2" "LSC" "Sludge" "Activated Sludge"], 'results_analyzer.xlsx', 'Sheet','sol.G','Range','A17:D17');
writematrix(["Preliminary 1"; "Preliminary 2"; "Preliminary 3"; "Primary 1"; "Primary 2"; "Primary 3"; "Secondary 1"; "Secondary 2"; "Secondary 3"; "Tertiary 1"; "Tertiary 2"; "Tertiary 3" ], 'results_analyzer.xlsx', 'Sheet','sol.G','Range','A18:A29');
writematrix([fG(5,1,2) fG(5,2,2) fG(5,3,2)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B18:D18');
writematrix([fG(6,1,2) fG(6,2,2) fG(6,3,2)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B19:D19');
writematrix([fG(7,1,2) fG(7,2,2) fG(7,3,2)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B20:D20');
writematrix([fG(9,1,2) fG(9,2,2) fG(9,3,2)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B21:D21');
writematrix([fG(10,1,2) fG(10,2,2) fG(10,3,2)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B22:D22');
writematrix([fG(11,1,2) fG(11,2,2) fG(11,3,2)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B23:D23');
writematrix([fG(13,1,2) fG(13,2,2) fG(13,3,2)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B24:D24');
writematrix([fG(14,1,2) fG(14,2,2) fG(14,3,2)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B25:D25');
writematrix([fG(15,1,2) fG(15,2,2) fG(15,3,2)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B26:D26');
writematrix([fG(17,1,2) fG(17,2,2) fG(17,3,2)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B27:D27');
writematrix([fG(18,1,2) fG(18,2,2) fG(18,3,2)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B28:D28');
writematrix([fG(19,1,2) fG(19,2,2) fG(19,3,2)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B29:D29');

writematrix(["Period 3" "LSC" "Sludge" "Activated Sludge"], 'results_analyzer.xlsx', 'Sheet','sol.G','Range','A33:D33');
writematrix(["Preliminary 1"; "Preliminary 2"; "Preliminary 3"; "Primary 1"; "Primary 2"; "Primary 3"; "Secondary 1"; "Secondary 2"; "Secondary 3"; "Tertiary 1"; "Tertiary 2"; "Tertiary 3" ], 'results_analyzer.xlsx', 'Sheet','sol.G','Range','A34:A45');
writematrix([fG(5,1,3) fG(5,2,3) fG(5,3,3)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B34:D34');
writematrix([fG(6,1,3) fG(6,2,3) fG(6,3,3)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B35:D35');
writematrix([fG(7,1,3) fG(7,2,3) fG(7,3,3)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B36:D36');
writematrix([fG(9,1,3) fG(9,2,3) fG(9,3,3)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B37:D37');
writematrix([fG(10,1,3) fG(10,2,3) fG(10,3,3)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B38:D38');
writematrix([fG(11,1,3) fG(11,2,3) fG(11,3,3)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B39:D39');
writematrix([fG(13,1,3) fG(13,2,3) fG(13,3,3)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B40:D40');
writematrix([fG(14,1,3) fG(14,2,3) fG(14,3,3)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B41:D41');
writematrix([fG(15,1,3) fG(15,2,3) fG(15,3,3)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B42:D42');
writematrix([fG(17,1,3) fG(17,2,3) fG(17,3,3)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B43:D43');
writematrix([fG(18,1,3) fG(18,2,3) fG(18,3,3)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B44:D44');
writematrix([fG(19,1,3) fG(19,2,3) fG(19,3,3)],'results_analyzer.xlsx', 'Sheet','sol.G','Range','B45:D45');


%sol.H
writematrix(["1" "2" "3"], 'results_analyzer.xlsx', 'Sheet','sol.H','Range','B2:D2');
writematrix(["Preliminary 1"; "Preliminary 2"; "Preliminary 3"; "Primary 1"; "Primary 2"; "Primary 3"; "Secondary 1"; "Secondary 2"; "Secondary 3"; "Tertiary 1"; "Tertiary 2"; "Tertiary 3" ], 'results_analyzer.xlsx', 'Sheet','sol.H','Range','A3:A14');
writematrix([fH(5,1,1) fH(5,2,1) fH(5,3,1)],'results_analyzer.xlsx', 'Sheet','sol.H','Range','B3:D3');
writematrix([fH(6,1,1) fH(6,2,1) fH(6,3,1)],'results_analyzer.xlsx', 'Sheet','sol.H','Range','B4:D4');
writematrix([fH(7,1,1) fH(7,2,1) fH(7,3,1)],'results_analyzer.xlsx', 'Sheet','sol.H','Range','B5:D5');
writematrix([fH(9,1,1) fH(9,2,1) fH(9,3,1)],'results_analyzer.xlsx', 'Sheet','sol.H','Range','B6:D6');
writematrix([fH(10,1,1) fH(10,2,1) fH(10,3,1)],'results_analyzer.xlsx', 'Sheet','sol.H','Range','B7:D7');
writematrix([fH(11,1,1) fH(11,2,1) fH(11,3,1)],'results_analyzer.xlsx', 'Sheet','sol.H','Range','B8:D8');
writematrix([fH(13,1,1) fH(13,2,1) fH(13,3,1)],'results_analyzer.xlsx', 'Sheet','sol.H','Range','B9:D9');
writematrix([fH(14,1,1) fH(14,2,1) fH(14,3,1)],'results_analyzer.xlsx', 'Sheet','sol.H','Range','B10:D10');
writematrix([fH(15,1,1) fH(15,2,1) fH(15,3,1)],'results_analyzer.xlsx', 'Sheet','sol.H','Range','B11:D11');
writematrix([fH(17,1,1) fH(17,2,1) fH(17,3,1)],'results_analyzer.xlsx', 'Sheet','sol.H','Range','B12:D12');
writematrix([fH(18,1,1) fH(18,2,1) fH(18,3,1)],'results_analyzer.xlsx', 'Sheet','sol.H','Range','B13:D13');
writematrix([fH(19,1,1) fH(19,2,1) fH(19,3,1)],'results_analyzer.xlsx', 'Sheet','sol.H','Range','B14:D14');

%sol.J
writematrix(["Period 1" "LSC" "Sludge" "Activated Sludge"], 'results_analyzer.xlsx', 'Sheet','sol.J','Range','A2:D2');
writematrix(["Preliminary 1"; "Preliminary 2"; "Preliminary 3"; "Primary 1"; "Primary 2"; "Primary 3"; "Secondary 1"; "Secondary 2"; "Secondary 3"; "Tertiary 1"; "Tertiary 2"; "Tertiary 3" ], 'results_analyzer.xlsx', 'Sheet','sol.J','Range','A3:A14');
writematrix([fJ(5,1,1) fJ(5,2,1) fJ(5,3,1)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B3:D3');
writematrix([fJ(6,1,1) fJ(6,2,1) fJ(6,3,1)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B4:D4');
writematrix([fJ(7,1,1) fJ(7,2,1) fJ(7,3,1)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B5:D5');
writematrix([fJ(9,1,1) fJ(9,2,1) fJ(9,3,1)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B6:D6');
writematrix([fJ(10,1,1) fJ(10,2,1) fJ(10,3,1)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B7:D7');
writematrix([fJ(11,1,1) fJ(11,2,1) fJ(11,3,1)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B8:D8');
writematrix([fJ(13,1,1) fJ(13,2,1) fJ(13,3,1)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B9:D9');
writematrix([fJ(14,1,1) fJ(14,2,1) fJ(14,3,1)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B10:D10');
writematrix([fJ(15,1,1) fJ(15,2,1) fJ(15,3,1)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B11:D11');
writematrix([fJ(17,1,1) fJ(17,2,1) fJ(17,3,1)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B12:D12');
writematrix([fJ(18,1,1) fJ(18,2,1) fJ(18,3,1)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B13:D13');
writematrix([fJ(19,1,1) fJ(19,2,1) fJ(19,3,1)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B14:D14');

writematrix(["Period 2" "LSC" "Sludge" "Activated Sludge"], 'results_analyzer.xlsx', 'Sheet','sol.J','Range','A18:D18');
writematrix(["Preliminary 1"; "Preliminary 2"; "Preliminary 3"; "Primary 1"; "Primary 2"; "Primary 3"; "Secondary 1"; "Secondary 2"; "Secondary 3"; "Tertiary 1"; "Tertiary 2"; "Tertiary 3" ], 'results_analyzer.xlsx', 'Sheet','sol.J','Range','A19:A30');
writematrix([fJ(5,1,2) fJ(5,2,2) fJ(5,3,2)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B19:D19');
writematrix([fJ(6,1,2) fJ(6,2,2) fJ(6,3,2)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B20:D20');
writematrix([fJ(7,1,2) fJ(7,2,2) fJ(7,3,2)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B21:D21');
writematrix([fJ(9,1,2) fJ(9,2,2) fJ(9,3,2)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B22:D22');
writematrix([fJ(10,1,2) fJ(10,2,2) fJ(10,3,2)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B23:D23');
writematrix([fJ(11,1,2) fJ(11,2,2) fJ(11,3,2)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B24:D24');
writematrix([fJ(13,1,2) fJ(13,2,2) fJ(13,3,2)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B25:D25');
writematrix([fJ(14,1,2) fJ(14,2,2) fJ(14,3,2)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B26:D26');
writematrix([fJ(15,1,2) fJ(15,2,2) fJ(15,3,2)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B27:D27');
writematrix([fJ(17,1,2) fJ(17,2,2) fJ(17,3,2)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B28:D28');
writematrix([fJ(18,1,2) fJ(18,2,2) fJ(18,3,2)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B29:D29');
writematrix([fJ(19,1,2) fJ(19,2,2) fJ(19,3,2)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B30:D30');

writematrix(["Period 3" "LSC" "Sludge" "Activated Sludge"], 'results_analyzer.xlsx', 'Sheet','sol.J','Range','A34:D34');
writematrix(["Preliminary 1"; "Preliminary 2"; "Preliminary 3"; "Primary 1"; "Primary 2"; "Primary 3"; "Secondary 1"; "Secondary 2"; "Secondary 3"; "Tertiary 1"; "Tertiary 2"; "Tertiary 3" ], 'results_analyzer.xlsx', 'Sheet','sol.J','Range','A35:A46');
writematrix([fJ(5,1,3) fJ(5,2,3) fJ(5,3,3)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B35:D35');
writematrix([fJ(6,1,3) fJ(6,2,3) fJ(6,3,3)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B36:D36');
writematrix([fJ(7,1,3) fJ(7,2,3) fJ(7,3,3)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B37:D37');
writematrix([fJ(9,1,3) fJ(9,2,3) fJ(9,3,3)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B38:D38');
writematrix([fJ(10,1,3) fJ(10,2,3) fJ(10,3,3)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B39:D39');
writematrix([fJ(11,1,3) fJ(11,2,3) fJ(11,3,3)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B40:D40');
writematrix([fJ(13,1,3) fJ(13,2,3) fJ(13,3,3)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B41:D41');
writematrix([fJ(14,1,3) fJ(14,2,3) fJ(14,3,3)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B42:D42');
writematrix([fJ(15,1,3) fJ(15,2,3) fJ(15,3,3)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B43:D43');
writematrix([fJ(17,1,3) fJ(17,2,3) fJ(17,3,3)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B44:D44');
writematrix([fJ(18,1,3) fJ(18,2,3) fJ(18,3,3)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B45:D45');
writematrix([fJ(19,1,3) fJ(19,2,3) fJ(19,3,3)],'results_analyzer.xlsx', 'Sheet','sol.J','Range','B46:D46');

%sol.N
writematrix(["Period 1" "BOD" "COD" "N" "P"], 'results_analyzer.xlsx', 'Sheet','sol.N','Range','A1:E1');
writematrix(["Storage 1"; "Storage 2"; "Storage 3"; "Storage 4"; "Output";], 'results_analyzer.xlsx', 'Sheet','sol.N','Range','A2:A6');
writematrix([fN(4,1,1) fN(4,2,1) fN(4,3,1) fN(4,4,1)],'results_analyzer.xlsx', 'Sheet','sol.N','Range','B2:E2');
writematrix([fN(8,1,1) fN(8,2,1) fN(8,3,1) fN(8,4,1)],'results_analyzer.xlsx', 'Sheet','sol.N','Range','B3:E3');
writematrix([fN(12,1,1) fN(12,2,1) fN(12,3,1) fN(12,4,1)],'results_analyzer.xlsx', 'Sheet','sol.N','Range','B4:E4');
writematrix([fN(16,1,1) fN(16,2,1) fN(16,3,1) fN(16,4,1)],'results_analyzer.xlsx', 'Sheet','sol.N','Range','B5:E5');
writematrix([fN(20,1,1) fN(20,2,1) fN(20,3,1) fN(20,4,1)],'results_analyzer.xlsx', 'Sheet','sol.N','Range','B6:E6');

writematrix(["Period 2" "BOD" "COD" "N" "P"], 'results_analyzer.xlsx', 'Sheet','sol.N','Range','A10:E10');
writematrix(["Storage 1"; "Storage 2"; "Storage 3"; "Storage 4"; "Output";], 'results_analyzer.xlsx', 'Sheet','sol.N','Range','A11:A15');
writematrix([fN(4,1,2) fN(4,2,2) fN(4,3,2) fN(4,4,2)],'results_analyzer.xlsx', 'Sheet','sol.N','Range','B11:E11');
writematrix([fN(8,1,2) fN(8,2,2) fN(8,3,2) fN(8,4,2)],'results_analyzer.xlsx', 'Sheet','sol.N','Range','B12:E12');
writematrix([fN(12,1,2) fN(12,2,2) fN(12,3,2) fN(12,4,2)],'results_analyzer.xlsx', 'Sheet','sol.N','Range','B13:E13');
writematrix([fN(16,1,2) fN(16,2,2) fN(16,3,2) fN(16,4,2)],'results_analyzer.xlsx', 'Sheet','sol.N','Range','B14:E14');
writematrix([fN(20,1,2) fN(20,2,2) fN(20,3,2) fN(20,4,2)],'results_analyzer.xlsx', 'Sheet','sol.N','Range','B15:E15');

writematrix(["Period 3" "BOD" "COD" "N" "P"], 'results_analyzer.xlsx', 'Sheet','sol.N','Range','A19:E19');
writematrix(["Storage 1"; "Storage 2"; "Storage 3"; "Storage 4"; "Output";], 'results_analyzer.xlsx', 'Sheet','sol.N','Range','A20:A24');
writematrix([fN(4,1,3) fN(4,2,3) fN(4,3,3) fN(4,4,3)],'results_analyzer.xlsx', 'Sheet','sol.N','Range','B20:E20');
writematrix([fN(8,1,3) fN(8,2,3) fN(8,3,3) fN(8,4,3)],'results_analyzer.xlsx', 'Sheet','sol.N','Range','B21:E21');
writematrix([fN(12,1,3) fN(12,2,3) fN(12,3,3) fN(12,4,3)],'results_analyzer.xlsx', 'Sheet','sol.N','Range','B22:E22');
writematrix([fN(16,1,3) fN(16,2,3) fN(16,3,3) fN(16,4,3)],'results_analyzer.xlsx', 'Sheet','sol.N','Range','B23:E23');
writematrix([fN(20,1,3) fN(20,2,3) fN(20,3,3) fN(20,4,3)],'results_analyzer.xlsx', 'Sheet','sol.N','Range','B24:E24');

%sol.O
writematrix(["Period 1" "LSC" "Sludge" "Activated Sludge"], 'results_analyzer.xlsx', 'Sheet','sol.O','Range','A1:E1');
writematrix(["Storage 1"; "Storage 2"; "Storage 3"; "Storage 4"; "Output";], 'results_analyzer.xlsx', 'Sheet','sol.O','Range','A2:A6');
writematrix([fO(4,1,1) fO(4,2,1) fO(4,3,1)],'results_analyzer.xlsx', 'Sheet','sol.O','Range','B2:E2');
writematrix([fO(8,1,1) fO(8,2,1) fO(8,3,1)],'results_analyzer.xlsx', 'Sheet','sol.O','Range','B3:E3');
writematrix([fO(12,1,1) fO(12,2,1) fO(12,3,1)],'results_analyzer.xlsx', 'Sheet','sol.O','Range','B4:E4');
writematrix([fO(16,1,1) fO(16,2,1) fO(16,3,1)],'results_analyzer.xlsx', 'Sheet','sol.O','Range','B5:E5');
writematrix([fO(20,1,1) fO(20,2,1) fO(20,3,1)],'results_analyzer.xlsx', 'Sheet','sol.O','Range','B6:E6');

writematrix(["Period 2" "LSC" "Sludge" "Activated Sludge"], 'results_analyzer.xlsx', 'Sheet','sol.O','Range','A10:E10');
writematrix(["Storage 1"; "Storage 2"; "Storage 3"; "Storage 4"; "Output";], 'results_analyzer.xlsx', 'Sheet','sol.O','Range','A11:A15');
writematrix([fO(4,1,2) fO(4,2,2) fO(4,3,2)],'results_analyzer.xlsx', 'Sheet','sol.O','Range','B11:E11');
writematrix([fO(8,1,2) fO(8,2,2) fO(8,3,2)],'results_analyzer.xlsx', 'Sheet','sol.O','Range','B12:E12');
writematrix([fO(12,1,2) fO(12,2,2) fO(12,3,2)],'results_analyzer.xlsx', 'Sheet','sol.O','Range','B13:E13');
writematrix([fO(16,1,2) fO(16,2,2) fO(16,3,2)],'results_analyzer.xlsx', 'Sheet','sol.O','Range','B14:E14');
writematrix([fO(20,1,2) fO(20,2,2) fO(20,3,2)],'results_analyzer.xlsx', 'Sheet','sol.O','Range','B15:E15');

writematrix(["Period 3" "LSC" "Sludge" "Activated Sludge"], 'results_analyzer.xlsx', 'Sheet','sol.O','Range','A19:E19');
writematrix(["Storage 1"; "Storage 2"; "Storage 3"; "Storage 4"; "Output";], 'results_analyzer.xlsx', 'Sheet','sol.O','Range','A20:A24');
writematrix([fO(4,1,3) fO(4,2,3) fO(4,3,3)],'results_analyzer.xlsx', 'Sheet','sol.O','Range','B20:E20');
writematrix([fO(8,1,3) fO(8,2,3) fO(8,3,3)],'results_analyzer.xlsx', 'Sheet','sol.O','Range','B21:E21');
writematrix([fO(12,1,3) fO(12,2,3) fO(12,3,3)],'results_analyzer.xlsx', 'Sheet','sol.O','Range','B22:E22');
writematrix([fO(16,1,3) fO(16,2,3) fO(16,3,3)],'results_analyzer.xlsx', 'Sheet','sol.O','Range','B23:E23');
writematrix([fO(20,1,3) fO(20,2,3) fO(20,3,3)],'results_analyzer.xlsx', 'Sheet','sol.O','Range','B24:E24');

%sol.B
writematrix(["Period 1" "1" "2" "3"], 'results_analyzer.xlsx', 'Sheet','sol.B','Range','A1:D1');
writematrix(["Biogas"; "Digestate";], 'results_analyzer.xlsx', 'Sheet','sol.B','Range','A2:A3');
writematrix([fB(1,1) fB(1,2) fB(1,3)],'results_analyzer.xlsx', 'Sheet','sol.B','Range','B2:D2');
writematrix([fB(2,1) fB(2,2) fB(2,3)],'results_analyzer.xlsx', 'Sheet','sol.B','Range','B3:D3');

%sol.A
writematrix(["Period 1"], 'results_analyzer.xlsx', 'Sheet','sol.A','Range','B1:B1');
writematrix(["Peiod 1"; "Period 2"; "Period 3"], 'results_analyzer.xlsx', 'Sheet','sol.A','Range','A2:A4');
writematrix([fA(1)],'results_analyzer.xlsx', 'Sheet','sol.A','Range','B2:B2');
writematrix([fA(2)],'results_analyzer.xlsx', 'Sheet','sol.A','Range','B3:B3');
writematrix([fA(3)],'results_analyzer.xlsx', 'Sheet','sol.A','Range','B4:B4');