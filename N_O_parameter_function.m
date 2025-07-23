function [sol, fval, exitflag, output] = N_O_parameter_function(N,O,Eq,Eb)
    

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
    % cm = cccc*2;
    cs = readmatrix("Parameters.xlsx",'Sheet','cs(s)','Range','B1:B4'); %capacity for storage s
    % cs = ccca*2;
    dw = readmatrix("Parameters.xlsx",'Sheet','dw(t)','Range','B2:B4'); %demand of water
    M = [100100; 112000; 400000]; %M = max(cm,[],"all");
    domJ = [zeros(4,3); 1 0 0; 1 0 0; 1 0 0; zeros(1,3); 0 1 0; 0 1 0; 0 1 0; zeros(1,3); 0 0 1; 0 0 1; 0 0 1; zeros(5,3)];
    domR = [0 1 0 0 1; 1 1 0 0 0; 1 1 0 0 0];
    
    %PARAMETERS EFFECTIVENESS
    L(:,:,1) = readmatrix("Parameters.xlsx",'Sheet','L(ilt)','Range','B2:E4');
    L(:,:,2) = readmatrix("Parameters.xlsx",'Sheet','L(ilt)','Range','H2:K4'); 
    L(:,:,3) = readmatrix("Parameters.xlsx",'Sheet','L(ilt)','Range','N2:Q4'); 
    K(:,:,1) = readmatrix("Parameters.xlsx",'Sheet','K(imt)','Range','B2:D4'); 
    K(:,:,2) = readmatrix("Parameters.xlsx",'Sheet','K(imt)','Range','G2:I4'); 
    K(:,:,3) = readmatrix("Parameters.xlsx",'Sheet','K(imt)','Range','L2:N4'); 
    % Eq = readmatrix("Parameters.xlsx",'Sheet','Eq(jl)','Range','B2:E21'); 
    % Eb = readmatrix("Parameters.xlsx",'Sheet','Eb(jm)','Range','B2:D21'); 
    P = readmatrix("Parameters.xlsx",'Sheet','P(l)','Range','B2:B5'); 
    W = readmatrix("Parameters.xlsx",'Sheet','W(m)','Range','B2:B4'); 
    
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
    % N = optimvar('N', j, l, t, 'LowerBound', 0, 'UpperBound', inf); %remaining percent quality parameter in the end of treatment
    % O = optimvar('O', j, m, t, 'LowerBound', 0, 'UpperBound', inf); %remaining percent byproducts in the end of treatment
    JJ = optimvar('JJ', m, t, 'LowerBound', 0, 'UpperBound', inf); %%dump J
    RR = optimvar('RR', m, n, t, 'LowerBound', 0, 'UpperBound', inf);
    
    %SYSTEM VARIABLES BDS
    B = optimvar('B', o, t, 'Type', 'continuous');
    A = optimvar('A', t, 'Type', 'continuous');
    udl = optimvar('udl', p, t, 'LowerBound', 0, 'UpperBound', inf);
    udc = optimvar('udc', q, t, 'LowerBound', 0, 'UpperBound', inf);
    udu = optimvar('udu', t, 'LowerBound', 0, 'UpperBound', inf);
    

    %OBJ VARIABLE
    ENVI = optimvar('ENVI', 'Type', 'continuous');

    %RANDOM VALUES
    randnum = readmatrix("Parameters.xlsx",'Sheet','RANDOM VALUES (For Monte Carlo)','Range','B1:B2000');
    
    WWTPNL = optimproblem;
    
    %%%%%%%%%%WWTP SUPPLY CHAIN CONSTRAINTS%%%%%%%%%%%%%%%%%%
    WWTPNL.Constraints.neta = e(1,4,:) == 1;
    WWTPNL.Constraints.netb = e(1,4,:) >= b(4,1,:);
    WWTPNL.Constraints.netc = b(4,1,:) >= e(2,8,:);
    WWTPNL.Constraints.netd = e(2,8,:) >= b(8,2,:) + b(8,3,:) + b(8,4,:);
    WWTPNL.Constraints.nete = b(8,2,:) >= e(3,12,:);
    WWTPNL.Constraints.netf = e(3,12,:) >= b(12,3,:) + b(12,4,:);
    WWTPNL.Constraints.netg = b(8,3,:) + b(12,3,:) >= e(4,16,:);
    WWTPNL.Constraints.neth = e(4,16,:) >= b(16,4,:);
    WWTPNL.Constraints.neti = b(8,4,:) + b(12,4,:) + b(16,4,:) >= e(5,20,:);
    WWTPNL.Constraints.netj = e(5,20,:) >= 1;
    
    %clustering (to treatments) constraints 
    WWTPNL.Constraints.tocluster1 = b(4,1,:) >= a(4,5,:);
    WWTPNL.Constraints.tocluster2 = b(4,1,:) >= a(4,6,:);
    WWTPNL.Constraints.tocluster3 = b(4,1,:) >= a(4,7,:);
    WWTPNL.Constraints.tocluster4 = b(8,2,:) >= a(8,9,:);
    WWTPNL.Constraints.tocluster5 = b(8,2,:) >= a(8,10,:);
    WWTPNL.Constraints.tocluster6 = b(8,2,:) >= a(8,11,:);
    WWTPNL.Constraints.tocluster7 = b(8,3,:) >= a(8,13,:);
    WWTPNL.Constraints.tocluster8 = b(8,3,:) >= a(8,14,:);
    WWTPNL.Constraints.tocluster9 = b(8,3,:) >= a(8,15,:);
    WWTPNL.Constraints.tocluster10 = b(8,4,:) >= a(8,17,:);
    WWTPNL.Constraints.tocluster11 = b(8,4,:) >= a(8,18,:);
    WWTPNL.Constraints.tocluster12 = b(8,4,:) >= a(8,19,:);
    WWTPNL.Constraints.tocluster13 = b(12,3,:) >= a(12,13,:);
    WWTPNL.Constraints.tocluster14 = b(12,3,:) >= a(12,14,:);
    WWTPNL.Constraints.tocluster15 = b(12,3,:) >= a(12,15,:);
    WWTPNL.Constraints.tocluster16 = b(12,4,:) >= a(12,17,:);
    WWTPNL.Constraints.tocluster17 = b(12,4,:) >= a(12,18,:);
    WWTPNL.Constraints.tocluster18 = b(12,4,:) >= a(12,19,:);
    WWTPNL.Constraints.tocluster19 = b(16,4,:) >= a(16,17,:);
    WWTPNL.Constraints.tocluster20 = b(16,4,:) >= a(16,18,:);
    WWTPNL.Constraints.tocluster21 = b(16,4,:) >= a(16,19,:);
    WWTPNL.Constraints.tocluster22 = sum(b(4,:,:),2) == 1;
    WWTPNL.Constraints.tocluster23 = sum(b(8,:,:),2) == 1;
    WWTPNL.Constraints.tocluster24 = sum(b(12,:,:),2) == 1;
    WWTPNL.Constraints.tocluster25 = sum(b(16,:,:),2) == 1;
    
    %clustering constraints (to storage)
    WWTPNL.Constraints.fromcluster1 = e(1,4,:) >= a(1,4,:);
    WWTPNL.Constraints.fromcluster2 = e(1,4,:) >= a(2,4,:);
    WWTPNL.Constraints.fromcluster3 = e(1,4,:) >= a(3,4,:);
    WWTPNL.Constraints.fromcluster4 = e(2,8,:) >= a(5,8,:);
    WWTPNL.Constraints.fromcluster5 = e(2,8,:) >= a(6,8,:);
    WWTPNL.Constraints.fromcluster6 = e(2,8,:) >= a(7,8,:);
    WWTPNL.Constraints.fromcluster7 = e(3,12,:) >= a(9,12,:);
    WWTPNL.Constraints.fromcluster8 = e(3,12,:) >= a(10,12,:);
    WWTPNL.Constraints.fromcluster9 = e(3,12,:) >= a(11,12,:);
    WWTPNL.Constraints.fromcluster10 = e(4,16,:) >= a(13,16,:);
    WWTPNL.Constraints.fromcluster11 = e(4,16,:) >= a(14,16,:);
    WWTPNL.Constraints.fromcluster12 = e(4,16,:) >= a(15,16,:);
    WWTPNL.Constraints.fromcluster13 = e(5,20,:) >= a(17,20,:);
    WWTPNL.Constraints.fromcluster14 = e(5,20,:) >= a(18,20,:);
    WWTPNL.Constraints.fromcluster15 = e(5,20,:) >= a(19,20,:);
    WWTPNL.Constraints.fromcluster16 = sum(e(:,4,:),1) == 1;
    WWTPNL.Constraints.fromcluster17 = sum(e(:,8,:),1) == 1;
    WWTPNL.Constraints.fromcluster18 = sum(e(:,12,:),1) == 1;
    WWTPNL.Constraints.fromcluster19 = sum(e(:,16,:),1) == 1;
    WWTPNL.Constraints.fromcluster20 = sum(e(:,20,:),1) == 1;
    
    %NEW
    WWTPNL.Constraints.binactive1 = x(:,:,1) <= M(1)*a(:,:,1);
    WWTPNL.Constraints.binactive2 = x(:,:,2) <= M(2)*a(:,:,2);
    WWTPNL.Constraints.binactive3 = x(:,:,3) <= M(3)*a(:,:,3);
    %SP.Constraints.asdfasd = x <= M*a;
    
    %input
    WWTPNL.Constraints.input1 = f(1,:)' == squeeze(x(1,4,:));
    WWTPNL.Constraints.input2 = f(2,:)' == squeeze(x(2,4,:));
    WWTPNL.Constraints.input3 = f(3,:)' == squeeze(x(3,4,:));
      
    %flow
    WWTPNL.Constraints.flow1 = x(1,4,:) + x(2,4,:) + x(3,4,:) == x(4,5,:) + x(4,6,:) + x(4,7,:);
    WWTPNL.Constraints.flow2 = x(4,5,:) >= x(5,8,:);
    WWTPNL.Constraints.flow3 = x(4,6,:) >= x(6,8,:);
    WWTPNL.Constraints.flow4 = x(4,7,:) >= x(7,8,:);
    WWTPNL.Constraints.flow5 = x(5,8,:) + x(6,8,:) + x(7,8,:) == x(8,9,:) + x(8,10,:) + x(8,11,:) + x(8,13,:) + x(8,14,:) + x(8,15,:) + x(8,17,:) + x(8,18,:) + x(8,19,:);
    WWTPNL.Constraints.flow6 = x(8,9,:) >= x(9,12,:);
    WWTPNL.Constraints.flow7 = x(8,10,:) >= x(10,12,:);
    WWTPNL.Constraints.flow8 = x(8,11,:) >= x(11,12,:);
    WWTPNL.Constraints.flow9 = x(9,12,:) + x(10,12,:) + x(11,12,:) == x(12,13,:) + x(12,14,:) + x(12,15,:) + x(12,17,:) + x(12,18,:) + x(12,19,:);
    WWTPNL.Constraints.flow10 = x(8,13,:) + x(12,13,:) >= x(13,16,:);
    WWTPNL.Constraints.flow11 = x(8,14,:) + x(12,14,:) >= x(14,16,:);
    WWTPNL.Constraints.flow12 = x(8,15,:) + x(12,15,:) >= x(15,16,:);
    WWTPNL.Constraints.flow13 = x(13,16,:) + x(14,16,:) + x(15,16,:) == x(16,17,:) + x(16,18,:) + x(16,19,:);
    WWTPNL.Constraints.flow14 = x(8,17,:) + x(12,17,:) + x(16,17,:) >= x(17,20,:);
    WWTPNL.Constraints.flow15 = x(8,18,:) + x(12,18,:) + x(16,18,:) >= x(18,20,:);
    WWTPNL.Constraints.flow16 = x(8,19,:) + x(12,19,:) + x(16,19,:) >= x(19,20,:);
    
    %demand 
    % WWTPNL.Constraints.demand1 = squeeze(x(17,20,:) + x(18,20,:) + x(19,20,:)) >= dw;
    
    %machine capacity
    WWTPNL.Constraints.machinecapacity1 = squeeze(x(4,5,:)) <= repmat(cm(1),[t,1]);
    WWTPNL.Constraints.machinecapacity2 = squeeze(x(4,6,:)) <= repmat(cm(2),[t,1]);
    WWTPNL.Constraints.machinecapacity3 = squeeze(x(4,7,:)) <= repmat(cm(3),[t,1]);
    WWTPNL.Constraints.machinecapacity4 = squeeze(x(8,9,:)) <= repmat(cm(4),[t,1]);
    WWTPNL.Constraints.machinecapacity5 = squeeze(x(8,10,:)) <= repmat(cm(5),[t,1]);
    WWTPNL.Constraints.machinecapacity6 = squeeze(x(8,11,:)) <= repmat(cm(6),[t,1]);
    WWTPNL.Constraints.machinecapacity7 = squeeze(x(8,13,:) + x(12,13,:)) <= repmat(cm(7),[t,1]);
    WWTPNL.Constraints.machinecapacity8 = squeeze(x(8,14,:) + x(12,14,:)) <= repmat(cm(8),[t,1]);
    WWTPNL.Constraints.machinecapacity9 = squeeze(x(8,15,:) + x(12,15,:)) <= repmat(cm(9),[t,1]);
    WWTPNL.Constraints.machinecapacity10 = squeeze(x(8,17,:) + x(12,17,:) + x(16,17,:)) <= repmat(cm(10),[t,1]);
    WWTPNL.Constraints.machinecapacity11 = squeeze(x(8,18,:) + x(12,18,:) + x(16,18,:)) <= repmat(cm(11),[t,1]);
    WWTPNL.Constraints.machinecapacity12 = squeeze(x(8,19,:) + x(12,19,:) + x(16,19,:)) <= repmat(cm(12),[t,1]);
    
    %storage capacity
    WWTPNL.Constraints.storagecapacity1 = squeeze(x(1,4,:) + x(2,4,:) + x(3,4,:)) <= repmat(cs(1),[t,1]);
    WWTPNL.Constraints.storagecapacity2 = squeeze(x(5,8,:) + x(6,8,:) + x(7,8,:)) <= repmat(cs(2),[t,1]);
    WWTPNL.Constraints.storagecapacity3 = squeeze(x(9,12,:) + x(10,12,:) + x(11,12,:)) <= repmat(cs(3),[t,1]);
    WWTPNL.Constraints.storagecapacity4 = squeeze(x(13,16,:) + x(14,16,:) + x(15,16,:)) <= repmat(cs(4),[t,1]);
    
    %%%%%%%%%%BYPRODUCT DISPOSAL SYSTEM CONSTRAINTS%%%%%%%%%%%%%%%%%%
    WWTPNL.Constraints.BDS1 =  squeeze(sum(J,1)) == squeeze(sum(R,2));
    WWTPNL.Constraints.BDS2 = squeeze(sum(R(:,1,:),1)) .* repmat(Q(1),[1,t])' == B(1,:)';
    WWTPNL.Constraints.BDS3 = squeeze(sum(R(:,1,:),1)) .* repmat(Q(2),[1,t])' == B(2,:)';
    WWTPNL.Constraints.BDS4 = B >= squeeze(sum(S,2));
    WWTPNL.Constraints.BDS5 = squeeze(sum(R(:,2,:),1)) + squeeze(sum(S(:,2,:),1)) == sum(T,1)';
    WWTPNL.Constraints.BDS6 = squeeze(sum(R(:,3,:),1)) + squeeze(sum(S(:,3,:),1)) == sum(U,1)';
    WWTPNL.Constraints.BDS7 =     T + udl >= dl;
    WWTPNL.Constraints.BDS8 =     U + udc >= dc;
    WWTPNL.Constraints.BDS9 =     squeeze(sum(R,1)) + squeeze(sum(S,1)) <= repmat(cp,[1,t]);
    WWTPNL.Constraints.BDS10 =     squeeze(sum(R(:,5,:),1)) + squeeze(sum(S(:,5,:),1)) == sum(V,1)';
    WWTPNL.Constraints.BDS11 =     A + udu >= du;
    WWTPNL.Constraints.BDS12 =     (squeeze(sum(S(:,4,:),1)) .* repmat(C,[1,t])') + (sum(V .* repmat(D,[1,t]) ,1)') == A;
    
    %%%%%%%%%%QUALITY CONSTRAINTS - NON LINEAR%%%%%%%%%%%%%%%%%%
    %computation of average input percentage for qp (aipq)
    aipqnl(1,:) = sum(f .* squeeze(L(:,1,:)),1) ./ sum(f,1);
    aipqnl(2,:) = sum(f .* squeeze(L(:,2,:)),1) ./ sum(f,1);
    aipqnl(3,:) = sum(f .* squeeze(L(:,3,:)),1) ./ sum(f,1);
    aipqnl(4,:) = sum(f .* squeeze(L(:,4,:)),1) ./ sum(f,1);

    %computation of average input percentage for bp (aipb)
    aipbnl(1,:) = sum(f .* squeeze(K(:,1,:)),1) ./ sum(f,1);
    aipbnl(2,:) = sum(f .* squeeze(K(:,2,:)),1) ./ sum(f,1);
    aipbnl(3,:) = sum(f .* squeeze(K(:,3,:)),1) ./ sum(f,1);

    % WWTPNL.Constraints.askdhf = aipqnl == squeeze(N(4,:,:));
    % WWTPNL.Constraints.askdha = aipbnl == squeeze(O(4,:,:));

    %-------FROM STORAGE 1--------%
    prelimqpvol1 = optimconstr(j,l,t);
    for j = 5:7
        for l = 1:l
            prelimqpvol1(j,l,:) = squeeze(x(4,j,:)) .* aipqnl(l,:)' .* (1 - squeeze(repmat(Eq(j,l),[1,1,t]))) == squeeze(F(j,l,:));
        end
    end
    WWTPNL.Constraints.fromstorage1 = prelimqpvol1;

    prelimbpvol1 = optimconstr(j,m,t);
    for j = 5:7
        for m = 1:m
            prelimbpvol1(j,m,:) = squeeze(x(4,j,:)) .* aipbnl(m,:)' .* (1 - squeeze(repmat(Eb(j,m),[1,1,t]))) == squeeze(G(j,m,:));
        end
    end
    WWTPNL.Constraints.fromstorage3 = prelimbpvol1;

    prelimbpvol1tobds = optimconstr(j,m,t);
    for j = 5:7
        for m = 1:m
            prelimbpvol1tobds(j,m,:) = squeeze(x(4,j,:)) .* aipbnl(m,:)' .* (squeeze(repmat(Eb(j,m),[1,1,t]))) == squeeze(J(j,m,:));
        end
    end
    WWTPNL.Constraints.fromstorage3aa = prelimbpvol1tobds;

    prelimclean1 = optimconstr(j,t);
    stor2total1 = optimconstr(j,t);
    for j = 5:7
        prelimclean1(j,:) = squeeze(x(4,j,:)) .* (1 - sum(aipqnl,1)' - sum(aipbnl,1)') == H(j,:)';
        stor2total1(j,:) = squeeze(x(j,8,:)) == squeeze(sum(F(j,:,:),2)) + squeeze(sum(G(j,:,:),2)) + H(j,:)';
    end
    WWTPNL.Constraints.fromstorage5 = prelimclean1;
    WWTPNL.Constraints.fromstorage7 = stor2total1;

    %-------FROM STORAGE 2--------%
    fromstorage2qpvol1 = optimconstr(j,l,t);
    for l = 1:l
        for j = 9:11
            fromstorage2qpvol1(j,l,:) = squeeze(x(8,j,:)) .* squeeze(N(8,l,:)) .* (1 - squeeze(repmat(Eq(j,l),[1,1,t]))) == squeeze(F(j,l,:));
        end
        for j = 13:15
            fromstorage2qpvol1(j,l,:) = (squeeze(x(8,j,:)) .* squeeze(N(8,l,:)) .* (1 - squeeze(repmat(Eq(j,l),[1,1,t])))) + (squeeze(x(12,j,:)) .* squeeze(N(12,l,:)) .* (1 - squeeze(repmat(Eq(j,l),[1,1,t])))) == squeeze(F(j,l,:));
        end
        for j = 17:19
            fromstorage2qpvol1(j,l,:) = (squeeze(x(8,j,:)) .* squeeze(N(8,l,:)) .* (1 - squeeze(repmat(Eq(j,l),[1,1,t])))) + (squeeze(x(12,j,:)) .* squeeze(N(12,l,:)) .* (1 - squeeze(repmat(Eq(j,l),[1,1,t])))) + (squeeze(x(16,j,:)) .* squeeze(N(16,l,:)) .* (1 - squeeze(repmat(Eq(j,l),[1,1,t])))) == squeeze(F(j,l,:));
        end
    end
    WWTPNL.Constraints.fromstorage13 = fromstorage2qpvol1;

    fromstorage2bpvol1 = optimconstr(j,m,t);
    for m = 1:m
        for j = 9:11
            fromstorage2bpvol1(j,m,:) = squeeze(x(8,j,:)) .* squeeze(O(8,m,:)) .* (1 - squeeze(repmat(Eb(j,m),[1,1,t]))) == squeeze(G(j,m,:));
        end
        for j = 13:15
            fromstorage2bpvol1(j,m,:) = (squeeze(x(8,j,:)) .* squeeze(O(8,m,:)) .* (1 - squeeze(repmat(Eb(j,m),[1,1,t])))) + (squeeze(x(12,j,:)) .* squeeze(O(12,m,:)) .* (1 - squeeze(repmat(Eb(j,m),[1,1,t])))) == squeeze(G(j,m,:));
        end
        for j = 17:19
            fromstorage2bpvol1(j,m,:) = (squeeze(x(8,j,:)) .* squeeze(O(8,m,:)) .* (1 - squeeze(repmat(Eb(j,m),[1,1,t])))) + (squeeze(x(12,j,:)) .* squeeze(O(12,m,:)) .* (1 - squeeze(repmat(Eb(j,m),[1,1,t])))) + (squeeze(x(16,j,:)) .* squeeze(O(16,m,:)) .* (1 - squeeze(repmat(Eb(j,m),[1,1,t])))) == squeeze(G(j,m,:));
        end
    end
    WWTPNL.Constraints.fromstorage15 = fromstorage2bpvol1;

    fromstorage2bdsvol1 = optimconstr(j,m,t);
    for m = 1:m
        for j = 9:11
            fromstorage2bdsvol1(j,m,:) = squeeze(x(8,j,:)) .* squeeze(O(8,m,:)) .* squeeze(repmat(Eb(j,m),[1,1,t])) == squeeze(J(j,m,:));
        end
        for j = 13:15
            fromstorage2bdsvol1(j,m,:) = (squeeze(x(8,j,:)) .* squeeze(O(8,m,:)) .* squeeze(repmat(Eb(j,m),[1,1,t]))) + (squeeze(x(12,j,:)) .* squeeze(O(12,m,:)) .* squeeze(repmat(Eb(j,m),[1,1,t]))) == squeeze(J(j,m,:));
        end
        for j = 17:19
            fromstorage2bdsvol1(j,m,:) = (squeeze(x(8,j,:)) .* squeeze(O(8,m,:)) .* squeeze(repmat(Eb(j,m),[1,1,t]))) + (squeeze(x(12,j,:)) .* squeeze(O(12,m,:)) .* squeeze(repmat(Eb(j,m),[1,1,t]))) + (squeeze(x(16,j,:)) .* squeeze(O(16,m,:)) .* squeeze(repmat(Eb(j,m),[1,1,t]))) == squeeze(J(j,m,:));
        end
    end
    WWTPNL.Constraints.fromstorage17 = fromstorage2bdsvol1;
    
    fromstor2clean1 = optimconstr(j,t);
    for j = 9:11
        fromstor2clean1(j,:) = squeeze(x(8,j,:)) .* (1 - squeeze(sum(N(8,:,:),2)) - squeeze(sum(O(8,:,:),2))) == H(j,:)';
    end
    for j = 13:15
        fromstor2clean1(j,:) = (squeeze(x(8,j,:)) .* (1 - squeeze(sum(N(8,:,:),2)) - squeeze(sum(O(8,:,:),2)))) + (squeeze(x(12,j,:)) .* (1 - squeeze(sum(N(12,:,:),2)) - squeeze(sum(O(12,:,:),2)))) == H(j,:)';
    end
    for j = 17:19
        fromstor2clean1(j,:) = (squeeze(x(8,j,:)) .* (1 - squeeze(sum(N(8,:,:),2)) - squeeze(sum(O(8,:,:),2)))) + (squeeze(x(12,j,:)) .* (1 - squeeze(sum(N(12,:,:),2)) - squeeze(sum(O(12,:,:),2)))) + (squeeze(x(16,j,:)) .* (1 - squeeze(sum(N(16,:,:),2)) - squeeze(sum(O(16,:,:),2)))) == H(j,:)';
    end
    WWTPNL.Constraints.fromstorage19 = fromstor2clean1;

    %total in storage 3
    stor3total1 = optimconstr(j,t);
    for j = 9:11
        stor3total1(j,:) = squeeze(x(j,12,:)) == squeeze(sum(F(j,:,:),2)) + squeeze(sum(G(j,:,:),2)) + H(j,:)';
    end
    WWTPNL.Constraints.fromstorage21 = stor3total1;

    %-------FROM STORAGE 3--------%
    %total in storage 4
    stor4total1 = optimconstr(j,t);
    for j = 13:15
        stor4total1(j,:) = squeeze(x(j,16,:)) == squeeze(sum(F(j,:,:),2)) + squeeze(sum(G(j,:,:),2)) + H(j,:)';
    end
    WWTPNL.Constraints.fromstorage27 = stor4total1;

    %-------FROM STORAGE 4--------%
    %total in output
    outputtotal1 = optimconstr(j,t);
    for j = 17:19
        outputtotal1(j,:) = squeeze(x(j,20,:)) == squeeze(sum(F(j,:,:),2)) + squeeze(sum(G(j,:,:),2)) + H(j,:)';
    end
    WWTPNL.Constraints.fromstorage33 = outputtotal1;

    %OBJECTIVE
    WWTPNL.Constraints.MULTOBJVAR = ENVI == sum(x .* repmat(PEw,[1,1,t]),'all') + sum(x .* repmat(TEw,[1,1,t]),'all')...
                                    + sum(R .* repmat(PEb,[1,1,t]),'all') + sum(R .* repmat(DEb,[1,1,t]),'all')...
                                    + sum(S .* repmat(PEn,[1,1,t]),'all') + sum(S .* repmat(DEn,[1,1,t]),'all')...
                                    + sum(V .* repmat(PEt,[1,t]),'all') + sum(U .* repmat(PEc,[1,t]),'all') + sum(T .* repmat(DEl,[1,t]),'all')...
                                    - sum(S .* repmat(GHG,[1,1,t]),'all');
    
    WWTPNL.ObjectiveSense = 'minimize';
    
    WWTPNL.Objective.COST = sum(x .* repmat(OCw,[1,1,t]),'all') + sum(x .* repmat(TCw,[1,1,t]),'all') + sum(a .* repmat(SC,[1,1,t]),'all') + sum(squeeze(sum(x,2)) .* repmat(HC,[1,t]),'all')... - sum(x .* repmat(SPw,[1,1,t]),'all')...
                         + sum(R .* repmat(OCb,[1,1,t]),'all') + sum(S .* repmat(OCn,[1,1,t]),'all') + sum(V .* OCt,'all')...
                         + sum(R .* repmat(TCb,[1,1,t]),'all') + sum(S .* repmat(TCn,[1,1,t]),'all') + sum(T .* TCl,'all') + sum(U .* TCc,'all')...
                         + sum(R .* repmat(DCb,[1,1,t]),'all') + sum(T .* DCl,'all') - sum(U .* SPc,'all') + ENVI*0.3; %+ sum(udu,'all')*10 + sum(udc,'all')*10 + sum(udl,'all')*10;
    
    [sol, fval, exitflag, output] = solve(WWTPNL);

end