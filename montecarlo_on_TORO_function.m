function [TOCW, xx, fval, MN, MO, TotCost, fENVI] = montecarlo_on_TORO_function(flow, Eq, Eb)
    
    %%%%%%%%%%%%%%%%RATIO CONVERSIONS%%%%%%%%%%%%%%%%
    flowratio = zeros(20,20,3);
    totprelim = flow(4,5,:) + flow(4,6,:) + flow(4,7,:);
    flowratio(4,5,:) = flow(4,5,:) ./ totprelim;
    flowratio(4,6,:) = flow(4,6,:) ./ totprelim;
    flowratio(4,7,:) = flow(4,7,:) ./ totprelim;
    
    totstor1toprim = flow(8,9,:) + flow(8,10,:) + flow(8,11,:);
    if totstor1toprim <= 0.1
        flowratio(8,9,:) = 0;
        flowratio(8,10,:) = 0;
        flowratio(8,11,:) = 0;
    else
        flowratio(8,9,:) = flow(8,9,:) ./ totstor1toprim;
        flowratio(8,10,:) = flow(8,10,:) ./ totstor1toprim;
        flowratio(8,11,:) = flow(8,11,:) ./ totstor1toprim;
    end
    
    totstor1tosec = flow(8,13,:) + flow(8,14,:) + flow(8,15,:);
    if totstor1tosec <= 0.1
        flowratio(8,13,:) = 0;
        flowratio(8,14,:) = 0;
        flowratio(8,15,:) = 0;
    else
        flowratio(8,13,:) = flow(8,13,:) ./ totstor1tosec;
        flowratio(8,14,:) = flow(8,14,:) ./ totstor1tosec;
        flowratio(8,15,:) = flow(8,15,:) ./ totstor1tosec;
    end
    
    totstor1toter = flow(8,17,:) + flow(8,18,:) + flow(8,19,:);
    if totstor1toter <= 0.1
        flowratio(8,17,:) = 0;
        flowratio(8,18,:) = 0;
        flowratio(8,19,:) = 0;
    else
        flowratio(8,17,:) = flow(8,17,:) ./ totstor1toter;
        flowratio(8,18,:) = flow(8,18,:) ./ totstor1toter;
        flowratio(8,19,:) = flow(8,19,:) ./ totstor1toter;
    end
    
    totstor2tosec = flow(12,13,:) + flow(12,14,:) + flow(12,15,:);
    if totstor2tosec <= 0.1
        flowratio(12,13,:) = 0;
        flowratio(12,14,:) = 0;
        flowratio(12,15,:) = 0;
    else
        flowratio(12,13,:) = flow(12,13,:) ./ totstor2tosec;
        flowratio(12,14,:) = flow(12,14,:) ./ totstor2tosec;
        flowratio(12,15,:) = flow(12,15,:) ./ totstor2tosec;
    end
    
    totstor2toter = flow(12,17,:) + flow(12,18,:) + flow(12,19,:);
    if totstor2toter <= 0.1
        flowratio(12,17,:) = 0;
        flowratio(12,18,:) = 0;
        flowratio(12,19,:) = 0;
    else
        flowratio(12,17,:) = flow(12,17,:) ./ totstor2toter;
        flowratio(12,18,:) = flow(12,18,:) ./ totstor2toter;
        flowratio(12,19,:) = flow(12,19,:) ./ totstor2toter;
    end
    
    totstor3toter = flow(16,17,:) + flow(16,18,:) + flow(16,19,:);
    if totstor3toter <= 0.1
        flowratio(16,17,:) = 0;
        flowratio(16,18,:) = 0;
        flowratio(16,19,:) = 0;
    else
        flowratio(16,17,:) = flow(16,17,:) ./ totstor3toter;
        flowratio(16,18,:) = flow(16,18,:) ./ totstor3toter;
        flowratio(16,19,:) = flow(16,19,:) ./ totstor3toter;
    end
    
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
    % Eq = readmatrix("Parameters.xlsx",'Sheet','Eq(jl)','Range','B2:E21'); 
    % Eb = readmatrix("Parameters.xlsx",'Sheet','Eb(jm)','Range','B2:D21'); 
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
    
    MONTERUN = optimproblem;

    %NEW
    MONTERUN.Constraints.binactive1 = x(:,:,1) <= M(1)*a(:,:,1);
    MONTERUN.Constraints.binactive2 = x(:,:,2) <= M(2)*a(:,:,2);
    MONTERUN.Constraints.binactive3 = x(:,:,3) <= M(3)*a(:,:,3);
 
    % %Byproduct disposal system
    MONTERUN.Constraints.BDS1 =    squeeze(sum(J,1)) == squeeze(sum(R,2));
    MONTERUN.Constraints.BDS2 =    squeeze(sum(R(:,1,:),1)) .* repmat(Q(1),[1,t])' == B(1,:)';
    MONTERUN.Constraints.BDS3 =    squeeze(sum(R(:,1,:),1)) .* repmat(Q(2),[1,t])' == B(2,:)';
    MONTERUN.Constraints.BDS4 =     B == squeeze(sum(S,2));
    MONTERUN.Constraints.BDS5 =     squeeze(sum(R(:,2,:),1)) + squeeze(sum(S(:,2,:),1)) == sum(T,1)';
    MONTERUN.Constraints.BDS6 =     squeeze(sum(R(:,3,:),1)) + squeeze(sum(S(:,3,:),1)) == sum(U,1)';
    MONTERUN.Constraints.BDS7 =     T + udl >= dl; 
    MONTERUN.Constraints.BDS8 =     U + udc >= dc;
    MONTERUN.Constraints.BDS9 =     squeeze(sum(R,1)) + squeeze(sum(S,1)) <= repmat(cp,[1,t]);
    MONTERUN.Constraints.BDS10 =     squeeze(sum(R(:,5,:),1)) + squeeze(sum(S(:,5,:),1)) == sum(V,1)';
    MONTERUN.Constraints.BDS11 =     A + udu >= du;
    MONTERUN.Constraints.BDS12 =     (squeeze(sum(S(:,4,:),1)) .* repmat(C,[1,t])') + (sum(V .* repmat(D,[1,t]) ,1)') == A;
    
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
    MONTERUN.Constraints.ina = f(1,:)' == squeeze(x(1,4,:));
    MONTERUN.Constraints.inb = f(2,:)' == squeeze(x(2,4,:));
    MONTERUN.Constraints.inc = f(3,:)' == squeeze(x(3,4,:));
    MONTERUN.Constraints.ind = (x(1,4,:) + x(2,4,:) + x(3,4,:)) .* flowratio(4,5,:) == x(4,5,:);
    MONTERUN.Constraints.ine = (x(1,4,:) + x(2,4,:) + x(3,4,:)) .* flowratio(4,6,:) == x(4,6,:);
    MONTERUN.Constraints.inf = (x(1,4,:) + x(2,4,:) + x(3,4,:)) .* flowratio(4,7,:) == x(4,7,:);
    
    %-------FROM STORAGE 1--------%
    prelimqpvol = optimconstr(j,l,t);
    for j = 5:7
        for l = 1:l
            prelimqpvol(j,l,:) = squeeze(x(4,j,:)) .* aipq(l,:)' .* (1 - squeeze(repmat(Eq(j,l),[1,1,t]))) == squeeze(F(j,l,:));
        end
    end
    MONTERUN.Constraints.fromstorage1 = prelimqpvol;
    
    prelimbpvol = optimconstr(j,m,t);
    for j = 5:7
        for m = 1:m
            prelimbpvol(j,m,:) = squeeze(x(4,j,:)) .* aipb(m,:)' .* (1 - squeeze(repmat(Eb(j,m),[1,1,t]))) == squeeze(G(j,m,:));
        end
    end
    MONTERUN.Constraints.fromstorage2 = prelimbpvol;
    
    prelimtobds = optimconstr(j,m,t);
    for j = 5:7
        for m = 1:m
            prelimtobds(j,m,:) = squeeze(x(4,j,:)) .* aipb(m,:)' .* (squeeze(repmat(Eb(j,m),[1,1,t]))) == squeeze(J(j,m,:));
        end
    end
    MONTERUN.Constraints.fromstorage3 = prelimtobds;
    
    prelimclean = optimconstr(j,t);
    for j = 5:7
        prelimclean(j,:) = squeeze(x(4,j,:)) .* (1 - sum(aipq,1)' - sum(aipb,1)') == H(j,:)';
    end
    MONTERUN.Constraints.fromstorage4 = prelimclean;
    
    stor2total = optimconstr(j,t);
    for j = 5:7
        stor2total(j,:) = squeeze(x(j,8,:)) == squeeze(sum(F(j,:,:),2)) + squeeze(sum(G(j,:,:),2)) + H(j,:)';
    end
    MONTERUN.Constraints.fromstorage5 = stor2total;
    
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
    MONTERUN.Constraints.fromstorage6 = frstor2dec;
    
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
    MONTERUN.Constraints.fromstorage7 = frstor2qpvol;
    
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
    MONTERUN.Constraints.fromstorage8 = frstor2bpvol;
    
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
    MONTERUN.Constraints.fromstorage9 = frstor2tobds;
    
    primclean = optimconstr(j,t);
    for j = 9:11
        primclean(j,:) = squeeze(x(8,j,:)) - squeeze(flowratio(8,j,:)) .* squeeze((sum(F(5,:,:),2) + sum(F(6,:,:),2) + sum(F(7,:,:),2) + sum(G(5,:,:),2) + sum(G(6,:,:),2) + sum(G(7,:,:),2))) == H(j,:)';
    end
    MONTERUN.Constraints.fromstorage10 = primclean;
    
    stor3total = optimconstr(j,t);
    for j = 9:11
        stor3total(j,:) = squeeze(x(j,12,:)) == squeeze(sum(F(j,:,:),2)) + squeeze(sum(G(j,:,:),2)) + H(j,:)';
    end
    MONTERUN.Constraints.fromstorage11 = stor3total;
    
    
    %-------FROM STORAGE 3--------%
    frstor3dec = optimconstr(j,j,t);
    for j = 13:15
        frstor3dec(12,j,:) = (x(9,12,:) + x(10,12,:) + x(11,12,:)) .* flowratio(12,j,:) == x(12,j,:);
    end
    for j = 17:19
        frstor3dec(12,j,:) = (x(9,12,:) + x(10,12,:) + x(11,12,:)) .* flowratio(12,j,:) == x(12,j,:);
    end
    MONTERUN.Constraints.fromstorage12 = frstor3dec;
    
    secclean = optimconstr(j,t);
    for j = 13:15
        secclean(j,:) = (squeeze(x(8,j,:)) - squeeze(flowratio(8,j,:)) .* squeeze((sum(F(5,:,:),2) + sum(F(6,:,:),2) + sum(F(7,:,:),2) + sum(G(5,:,:),2) + sum(G(6,:,:),2) + sum(G(7,:,:),2)))) + (squeeze(x(12,j,:)) - squeeze(flowratio(12,j,:)) .* squeeze((sum(F(9,:,:),2) + sum(F(10,:,:),2) + sum(F(11,:,:),2) + sum(G(9,:,:),2) + sum(G(10,:,:),2) + sum(G(11,:,:),2)))) == H(j,:)';
    end
    MONTERUN.Constraints.fromstorage13 = secclean;
    
    stor4total = optimconstr(j,t);
    for j = 13:15
        stor4total(j,:) = squeeze(x(j,16,:)) == squeeze(sum(F(j,:,:),2)) + squeeze(sum(G(j,:,:),2)) + H(j,:)';
    end
    MONTERUN.Constraints.fromstorage14 = stor4total;
    
    %-------FROM STORAGE 4--------%
    frstor4dec = optimconstr(j,j,t);
    for j = 17:19
        frstor4dec(16,j,:) = (x(13,16,:) + x(14,16,:) + x(15,16,:)) .* flowratio(16,j,:) == x(16,j,:);
    end
    MONTERUN.Constraints.fromstorage15 = frstor4dec;
    
    terclean = optimconstr(j,t);
    for j = 17:19
        terclean(j,:) = (squeeze(x(8,j,:)) - squeeze(flowratio(8,j,:)) .* squeeze((sum(F(5,:,:),2) + sum(F(6,:,:),2) + sum(F(7,:,:),2) + sum(G(5,:,:),2) + sum(G(6,:,:),2) + sum(G(7,:,:),2)))) + (squeeze(x(12,j,:)) - squeeze(flowratio(12,j,:)) .* squeeze((sum(F(9,:,:),2) + sum(F(10,:,:),2) + sum(F(11,:,:),2) + sum(G(9,:,:),2) + sum(G(10,:,:),2) + sum(G(11,:,:),2)))) + (squeeze(x(16,j,:)) - squeeze(flowratio(16,j,:)) .* squeeze((sum(F(13,:,:),2) + sum(F(14,:,:),2) + sum(F(15,:,:),2) + sum(G(13,:,:),2) + sum(G(14,:,:),2) + sum(G(15,:,:),2)))) == H(j,:)';
    end
    MONTERUN.Constraints.fromstorage16 = terclean;
    
    outputtotal = optimconstr(j,t);
    for j = 17:19
        outputtotal(j,:) = squeeze(x(j,20,:)) == squeeze(sum(F(j,:,:),2)) + squeeze(sum(G(j,:,:),2)) + H(j,:)';
    end
    MONTERUN.Constraints.fromstorage17 = outputtotal;
    
    MONTERUN.Constraints.MULTOBJVAR = ENVI == sum(x .* repmat(PEw,[1,1,t]),'all') + sum(x .* repmat(TEw,[1,1,t]),'all')...
                                        + sum(R .* repmat(PEb,[1,1,t]),'all') + sum(R .* repmat(DEb,[1,1,t]),'all')...
                                        + sum(S .* repmat(PEn,[1,1,t]),'all') + sum(S .* repmat(DEn,[1,1,t]),'all')...
                                        + sum(V .* repmat(PEt,[1,t]),'all') + sum(U .* repmat(PEc,[1,t]),'all') + sum(T .* repmat(DEl,[1,t]),'all')...
                                        - sum(S .* repmat(GHG,[1,1,t]),'all');
    
    MONTERUN.ObjectiveSense = 'minimize';
    
    MONTERUN.Objective.COST = sum(x .* repmat(OCw,[1,1,t]),'all') + sum(x .* repmat(TCw,[1,1,t]),'all') + sum(a .* repmat(SC,[1,1,t]),'all') + sum(squeeze(sum(x,2)) .* repmat(HC,[1,t]),'all')...
                             + sum(R .* repmat(OCb,[1,1,t]),'all') + sum(S .* repmat(OCn,[1,1,t]),'all') + sum(V .* OCt,'all')...
                             + sum(R .* repmat(TCb,[1,1,t]),'all') + sum(S .* repmat(TCn,[1,1,t]),'all') + sum(T .* TCl,'all') + sum(U .* TCc,'all')...
                             + sum(R .* repmat(DCb,[1,1,t]),'all') + sum(T .* DCl,'all') - sum(U .* SPc,'all')...
                             + ENVI*0.3;
    
    [sol, fval, exitflag, output] = solve(MONTERUN);
    
    TOCW = sum(sol.x .* repmat(OCw,[1,1,t]),'all');

    xx = sol.x;
    
    MN = zeros(20,4,3);
    for l = 1:l
        %node 8
        if sol.x(5,8,:) + sol.x(6,8,:) + sol.x(7,8,:) == 0
            MN(8,l,:) = 0;
        else
            MN(8,l,:) = (sol.F(5,l,:) + sol.F(6,l,:) + sol.F(7,l,:) ) ./ (sol.x(5,8,:) + sol.x(6,8,:) + sol.x(7,8,:));
        end
        %node 12
        if sol.x(9,12,:) + sol.x(10,12,:) + sol.x(11,12,:) == 0
            MN(12,l,:) = 0;
        else
            MN(12,l,:) = (sol.F(9,l,:) + sol.F(10,l,:) + sol.F(11,l,:) ) ./ (sol.x(9,12,:) + sol.x(10,12,:) + sol.x(11,12,:));
        end
        %node 16
        if sol.x(13,16,:) + sol.x(14,16,:) + sol.x(15,16,:) == 0
            MN(16,l,:) = 0;
        else
            MN(16,l,:) = (sol.F(13,l,:) + sol.F(14,l,:) + sol.F(15,l,:) ) ./ (sol.x(13,16,:) + sol.x(14,16,:) + sol.x(15,16,:));
        end
        %node 20
        if sol.x(17,20,:) + sol.x(18,20,:) + sol.x(19,20,:) == 0
            MN(20,l,:) = 0;
        else
            MN(20,l,:) = (sol.F(17,l,:) + sol.F(18,l,:) + sol.F(19,l,:) ) ./ (sol.x(17,20,:) + sol.x(18,20,:) + sol.x(19,20,:));
        end
    end
    
    MO = zeros(20,3,3);
    for m = 1:m
        %node 8
        if sol.x(5,8,:) + sol.x(6,8,:) + sol.x(7,8,:) == 0
            MO(8,m,:) = 0;
        else
            MO(8,m,:) = (sol.G(5,m,:) + sol.G(6,m,:) + sol.G(7,m,:) ) ./ (sol.x(5,8,:) + sol.x(6,8,:) + sol.x(7,8,:));
        end
        %node 12
        if sol.x(9,12,:) + sol.x(10,12,:) + sol.x(11,12,:) == 0
            MO(12,m,:) = 0;
        else
            MO(12,m,:) = (sol.G(9,m,:) + sol.G(10,m,:) + sol.G(11,m,:) ) ./ (sol.x(9,12,:) + sol.x(10,12,:) + sol.x(11,12,:));
        end
        %node 16
        if sol.x(13,16,:) + sol.x(14,16,:) + sol.x(15,16,:) == 0
            MO(16,m,:) = 0;
        else
            MO(16,m,:) = (sol.G(13,m,:) + sol.G(14,m,:) + sol.G(15,m,:) ) ./ (sol.x(13,16,:) + sol.x(14,16,:) + sol.x(15,16,:));
        end
        %node 20
        if sol.x(17,20,:) + sol.x(18,20,:) + sol.x(19,20,:) == 0
            MO(20,m,:) = 0;
        else
            MO(20,m,:) = (sol.G(17,m,:) + sol.G(18,m,:) + sol.G(19,m,:) ) ./ (sol.x(17,20,:) + sol.x(18,20,:) + sol.x(19,20,:));
        end
    end

    TotCost = sum(sol.x .* repmat(OCw,[1,1,t]),'all') + sum(sol.x .* repmat(TCw,[1,1,t]),'all') + sum(sol.a .* repmat(SC,[1,1,t]),'all') + sum(squeeze(sum(sol.x,2)) .* repmat(HC,[1,t]),'all')... 
                             + sum(sol.R .* repmat(OCb,[1,1,t]),'all') + sum(sol.S .* repmat(OCn,[1,1,t]),'all') + sum(sol.V .* OCt,'all')...
                             + sum(sol.R .* repmat(TCb,[1,1,t]),'all') + sum(sol.S .* repmat(TCn,[1,1,t]),'all') + sum(sol.T .* TCl,'all') + sum(sol.U .* TCc,'all')...
                             + sum(sol.R .* repmat(DCb,[1,1,t]),'all') + sum(sol.T .* DCl,'all');
    fENVI = sol.ENVI;

    
end