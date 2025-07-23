function [sol, fval, exitflag, output] = x_parameter_function(x,Eq,Eb)

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
    % x = optimvar('x', j, j, t, 'LowerBound', 0, 'UpperBound', inf);
    
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
    
    WWTPNL = optimproblem;
    
    %quality and byproduct requirement
    WWTPNL.Constraints.qpreq = squeeze(N(20,:,:)) <= repmat(P,[1,t]);
    WWTPNL.Constraints.bpreq = squeeze(O(20,:,:)) <= repmat(W,[1,t]);

    % WWTPNL.Constraints.trial = squeeze(sum(J .* repmat(domJ,[1,1,t]),1)) == JJ;
    % WWTPNL.Constraints.trial1 = R .* repmat(domR,[1,1,t]) == RR;
    
    %%%%%%%%%%BYPRODUCT DISPOSAL SYSTEM CONSTRAINTS%%%%%%%%%%%%%%%%%%
    WWTPNL.Constraints.BDS1 = squeeze(sum(J,1)) == squeeze(sum(R,2));
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
    
    WWTPNL.Constraints.askdhf = aipqnl == squeeze(N(4,:,:));
    WWTPNL.Constraints.askdha = aipbnl == squeeze(O(4,:,:));
    
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
    % stor2total1 = optimconstr(j,t);
    for j = 5:7
        prelimclean1(j,:) = squeeze(x(4,j,:)) .* (1 - sum(aipqnl,1)' - sum(aipbnl,1)') == H(j,:)';
        % stor2total1(j,:) = squeeze(x(j,8,:)) == squeeze(sum(F(j,:,:),2)) + squeeze(sum(G(j,:,:),2)) + H(j,:)';
    end
    WWTPNL.Constraints.fromstorage5 = prelimclean1;
    % WWTPNL.Constraints.fromstorage7 = stor2total1;
    
    %non linear starting this point
    %quality and byproducts [percent] in storage 2
    stor2percentqp1 = optimconstr(l,t);
    for l = 1:l
        stor2percentqp1(l,:) = squeeze(x(5,8,:) + x(6,8,:) + x(7,8,:)) .* squeeze(N(8,l,:)) == squeeze(F(5,l,:) + F(6,l,:) + F(7,l,:));
    end
    WWTPNL.Constraints.fromstorage9 = stor2percentqp1;

    stor2percentbp1 = optimconstr(m,t);
    for m = 1:m
        stor2percentbp1(m,:) = squeeze(x(5,8,:) + x(6,8,:) + x(7,8,:)) .* squeeze(O(8,m,:)) == squeeze(G(5,m,:) + G(6,m,:) + G(7,m,:));
    end
    WWTPNL.Constraints.fromstorage11 = stor2percentbp1;

    %quality and by products [percent] in storage 3
    stor3percentqp1 = optimconstr(l,t);
    for t = 1:t
        for l = 1:l
            stor3percentqp1(l,t) = squeeze(x(9,12,t) + x(10,12,t) + x(11,12,t)) .* squeeze(N(12,l,t)) == squeeze(F(9,l,t) + F(10,l,t) + F(11,l,t));
        end
    end
    WWTPNL.Constraints.fromstorage23 = stor3percentqp1;

    stor3percentbp1 = optimconstr(m,t);
    for t = 1:t
        for m = 1:m
            stor3percentbp1(m,t) = squeeze(x(9,12,t) + x(10,12,t) + x(11,12,t)) .* squeeze(O(12,m,t)) == squeeze(G(9,m,t) + G(10,m,t) + G(11,m,t));
        end
    end
    WWTPNL.Constraints.fromstorage25 = stor3percentbp1;

    %quality and by products [percent] in storage 4
    stor4percentqp1 = optimconstr(l,t);
    for t = 1:t
        for l = 1:l
            stor4percentqp1(l,t) = squeeze(x(13,16,t) + x(14,16,t) + x(15,16,t)) .* squeeze(N(16,l,t)) == squeeze(F(13,l,t) + F(14,l,t) + F(15,l,t));
        end
    end
    WWTPNL.Constraints.fromstorage29 = stor4percentqp1;

    stor4percentbp1 = optimconstr(m,t);
    for t = 1:t
        for m = 1:m
            stor4percentbp1(m,t) = squeeze(x(13,16,t) + x(14,16,t) + x(15,16,t)) .* squeeze(O(16,m,t)) == squeeze(G(13,m,t) + G(14,m,t) + G(15,m,t));
        end
    end
    WWTPNL.Constraints.fromstorage31 = stor4percentbp1;

    %quality and by products [percent] in storage 4
    outputpercentqp1 = optimconstr(l,t);
    for l = 1:l
        outputpercentqp1(l,:) = squeeze(x(17,20,:) + x(18,20,:) + x(19,20,:)) .* squeeze(N(20,l,:)) == squeeze(F(17,l,:) + F(18,l,:) + F(19,l,:));
    end
    WWTPNL.Constraints.fromstorage35 = outputpercentqp1;

    outputpercentbp1 = optimconstr(m,t);
    for m = 1:m
        outputpercentbp1(m,:) = squeeze(x(17,20,:) + x(18,20,:) + x(19,20,:)) .* squeeze(O(20,m,:)) == squeeze(G(17,m,:) + G(18,m,:) + G(19,m,:));
    end
    WWTPNL.Constraints.fromstorage37 = outputpercentbp1;

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
                         + sum(R .* repmat(DCb,[1,1,t]),'all') + sum(T .* DCl,'all') - sum(U .* SPc,'all') + ENVI*0.3;% + sum(udu,'all')*10 + sum(udc,'all')*10 + sum(udl,'all')*10;
    
    [sol, fval, exitflag, output] = solve(WWTPNL);

end