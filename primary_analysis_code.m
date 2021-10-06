%% CLEAR WORKSPACE
clc;   
close all; 
clear;

%% GENERATE 3x3 DESIGNS

% Initialize values
sidenum = 3;
numDes = 200;
sel = 0.01; 
E = 1816200;
rnum = 250;
r = rnum*(10^-6);
CBasket = []; 
vfBasket = [];

% Generate array with nodal coordinates
NC = generateNC(sel,sidenum);

% Generate asymmetric 3x3 designs
asymmCA = gen_Des_func_V2(sidenum,sel,numDes,0.5,0,2);

% Generate symmetric 3x3 designs
symmCA = inputSymmCA();

% Generate flip symmetric 3x3 designs
vfsymmCA = gen_Des_func_V3(sidenum,sel,5,...
         0.5,0,2,2);
hfsymmCA = gen_Des_func_V3(sidenum,sel,5,...
         0.5,0,2,3);

% Populate 3x3 CABasket
CABasket_3x3 = [asymmCA';vfsymmCA;hfsymmCA;symmCA'];
for j = 1:1:length(CABasket_3x3)
    CA = cell2mat(CABasket_3x3(j));
    SortedCA = sortrows((sort(CA'))');
    CABasket_3x3(j) = {SortedCA};
end

% Isolate unique CAs
usedCAs = {}; uCAc = 0;
tobeusedCAs = {};
for v = 1:1:length(CABasket_3x3)
    % Initialize values
    CA = cell2mat(CABasket_3x3(v));
    
    % Check if current CA has been used already
    usedbool = 0;
    for k = 1:1:length(usedCAs)
        uCA = cell2mat(usedCAs(k));
        if size(CA,1) == size(uCA,1)
            if CA == uCA
                usedbool = 1;
                break
            end
        end
    end
    if usedbool == 0
        % Log current CA as used
        uCAc = uCAc + 1;
        usedCAs(uCAc) = {CA}; 
        tobeusedCAs(uCAc) = {CA};
    end
end
tobeusedCAs = tobeusedCAs(~cellfun('isempty',tobeusedCAs));

%% GENERATE 5x5 DESIGNS

% Initialize values
sidenum = 5;
numDes_as = 10;
numDes_fs = 10;
numDes_rs = 10;
sel = 0.01; 
E = 1816200;
rnum = 250;
r = rnum*(10^-6);
CBasket = []; 
vfBasket = [];

% Generate array with nodal coordinates
NC = generateNC(sel,sidenum);

% Generate asymmetric 5x5 designs
asymmCA_5x5 = gen_Des_func_V2(sidenum,sel,numDes_as,0.7,1,2);

% Restricting total count of all CA's considered
asymmCA_5x5 = asymmCA_5x5(1:600);
CABasket_3x3 = CABasket_3x3(1:600);

% Develop rotationally symmetric 5x5 designs 
rotsymmCA_5x5 = gen_Des_func_V3(sidenum,sel,numDes_rs,...
         0.7,0,2,1);

% Check rotationally symmetric 5x5 designs for repeatability
oldrotsymmCA = rotsymmCA_5x5;
rotsymmCA_5x5 = {}; goodcount = 1;
for i = 1:1:size(oldrotsymmCA,1)
    CA = cell2mat(oldrotsymmCA(i));
    repBool = repChecker_2D_V1(CA,sidenum);
    if repBool == 1
        rotsymmCA_5x5(goodcount) = {CA};
        goodcount = goodcount + 1;
    end
end

% Develop flip-symmetric 5x5 designs 
vflipsymmCA_5x5 = gen_Des_func_V3(sidenum,sel,numDes_fs,...
         0.7,0,2,2);
hflipsymmCA_5x5 = gen_Des_func_V3(sidenum,sel,numDes_fs,...
         0.7,0,2,3);    

% Populate full CABasket
CABasket_5x5 = [rotsymmCA_5x5';hflipsymmCA_5x5;vflipsymmCA_5x5;asymmCA_5x5];
for j = 1:1:length(CABasket_5x5)
    CA = cell2mat(CABasket_5x5(j));
    SortedCA = sortrows((sort(CA'))');
    CABasket_5x5(j) = {SortedCA};
end

% Isolate unique CAs
usedCAs = {}; uCAc = 0;
for v = 1:1:length(CABasket_5x5)
    % Initialize values
    CA = cell2mat(CABasket_5x5(v));
    
    % Check if current CA has been used already
    usedbool = 0;
    for k = 1:1:length(usedCAs)
        uCA = cell2mat(usedCAs(k));
        if size(CA,1) == size(uCA,1)
            if CA == uCA
                usedbool = 1;
                break
            end
        end
    end
    if usedbool == 0
        % Log current CA as used
        uCAc = uCAc + 1;
        usedCAs(uCAc) = {CA}; 
        tobeusedCAs(uCAc) = {CA};
    end
end
tobeusedCAs = tobeusedCAs(~cellfun('isempty',tobeusedCAs));


%% COMPUTE STIFFNESS VALUES
% Determine current grid size
if sidenum == 3
    CABasket = tobeusedCAs;
elseif sidenum == 5
    CABasket = tobeusedCAs;
end

% Loop through each design to solve for, plot C-matrix values
for cac = 1:1:size(CABasket,2)
    CA = cell2mat(CABasket(cac));
    
    % Calculate edge-related inputs
    [CA,edgelog] = edgeVals(CA,sidenum);
    
    % Calculate cross-sectional properties for edge members
    [csPropsA,csPropsB] = calcCSProps(r);
    
    % Determine nodal pairs for each constraint equation
    [NP,numLRCE,numTBCE] = findNodalPairs(sidenum,CA,NC);
    
    % Calculate C, volFrac
    [C,vf] = calcCTruss(NC,csPropsA,csPropsB,edgelog,NP,...
                             numLRCE,numTBCE,sidenum,sel,r,E,CA);

    % Record outputs for plotting/reference
    CBasket(:,:,cac) = C; 
    vfBasket = [vfBasket,vf];
    strout = strcat('CURRENT DESIGN: ',num2str(cac));
    disp(strout);
end

%% POSTPROCESS RESULTS
CBasket = CBasket.*100;
C11Basket = []; C12Basket = []; C21Basket = []; C22Basket = [];
C16Basket = []; C61Basket = []; C26Basket = []; C62Basket = [];
C66Basket = []; vxyBasket = []; vyxBasket = []; symmScores = [];
S11Basket = []; S22Basket = []; S12Basket = []; S21Basket = [];
sC11Basket = []; sC12Basket = []; sC21Basket = []; sC22Basket = [];
sC16Basket = []; sC61Basket = []; sC26Basket = []; sC62Basket = [];
sC66Basket = []; svxyBasket = []; svyxBasket = []; 
fsC11Basket = []; fsC12Basket = []; fsC21Basket = []; fsC22Basket = [];
fsC16Basket = []; fsC61Basket = []; fsC26Basket = []; fsC62Basket = [];
fsC66Basket = []; fsvxyBasket = []; fsvyxBasket = []; 
asymmL = 0; symmL = 0; fsymmL = 0; 
ascounter = 1; scounter = 1; fscounter = 1;
goodCABasket = {}; sgoodCABasket = {}; fsgoodCABasket = {}; 
skB = []; askB = []; fskB = [];
for k = 1:1:size(CBasket,3)
    Ctemp = CBasket(:,:,k); 
    Stemp = inv(Ctemp);
    vxy = -(Stemp(2,1)/Stemp(1,1));
    vyx = -(Stemp(1,2)/Stemp(2,2));
    miniC = Ctemp(1:2,1:2);
    Stemp2 = inv(miniC);
    vxy2 = -(Stemp2(2,1)/Stemp2(1,1));
    vyx2 = -(Stemp2(1,2)/Stemp2(2,2));
    
    % Calculate symmetry score
    tinyCA = cell2mat(CABasket(k));
    SS = symmHeuristic_2D(tinyCA,NC,sel);
    symmScores = [symmScores,SS]; 
    if SS == 0
        asymmL = asymmL + 1;
        vxyBasket = [vxyBasket,vxy]; 
        vyxBasket = [vyxBasket,vyx];
        C11Basket = [C11Basket,Ctemp(1,1)];
        C12Basket = [C12Basket,Ctemp(1,2)];
        C21Basket = [C21Basket,Ctemp(2,1)];
        C22Basket = [C22Basket,Ctemp(2,2)];
        C16Basket = [C16Basket,Ctemp(1,3)];
        C61Basket = [C61Basket,Ctemp(3,1)];
        C26Basket = [C26Basket,Ctemp(2,3)];
        C62Basket = [C62Basket,Ctemp(3,2)];
        C66Basket = [C66Basket,Ctemp(3,3)];
        goodCABasket(ascounter) = {tinyCA}; 
        ascounter = ascounter + 1;
        askB = [askB,k];
    elseif SS == 1
        fsymmL = fsymmL + 1;
        fsvxyBasket = [fsvxyBasket,vxy]; 
        fsvyxBasket = [fsvyxBasket,vyx];
        fsC11Basket = [fsC11Basket,Ctemp(1,1)];
        fsC12Basket = [fsC12Basket,Ctemp(1,2)];
        fsC21Basket = [fsC21Basket,Ctemp(2,1)];
        fsC22Basket = [fsC22Basket,Ctemp(2,2)];
        fsC16Basket = [fsC16Basket,Ctemp(1,3)];
        fsC61Basket = [fsC61Basket,Ctemp(3,1)];
        fsC26Basket = [fsC26Basket,Ctemp(2,3)];
        fsC62Basket = [fsC62Basket,Ctemp(3,2)];
        fsC66Basket = [fsC66Basket,Ctemp(3,3)];
        fsgoodCABasket(fscounter) = {tinyCA}; 
        fscounter = fscounter + 1;
        fskB = [fskB,k];
    elseif SS == 2
        symmL = symmL + 1;
        svxyBasket = [svxyBasket,vxy2]; 
        svyxBasket = [svyxBasket,vyx2];
        sC11Basket = [sC11Basket,Ctemp(1,1)];
        sC12Basket = [sC12Basket,Ctemp(1,2)];
        sC21Basket = [sC21Basket,Ctemp(2,1)];
        sC22Basket = [sC22Basket,Ctemp(2,2)];
        S11Basket = [S11Basket,Stemp(1,1)];
        S12Basket = [S12Basket,Stemp(1,2)];
        S21Basket = [S21Basket,Stemp(2,1)];
        S22Basket = [S22Basket,Stemp(2,2)];
        sC16Basket = [sC16Basket,Ctemp(1,3)];
        sC61Basket = [sC61Basket,Ctemp(3,1)];
        sC26Basket = [sC26Basket,Ctemp(2,3)];
        sC62Basket = [sC62Basket,Ctemp(3,2)];
        sC66Basket = [sC66Basket,Ctemp(3,3)];
        sgoodCABasket(scounter) = {tinyCA}; 
        scounter = scounter + 1;
        skB = [skB,k];
    end
end

%% PLOT RESULTS- 2D SCATTERPLOT (Shear-vF)
plot([vfBasket(skB),vfBasket(skB)],[sC66Basket,sC66Basket]./E,'Marker','.','MarkerSize',8,'LineStyle','none','Color',[0.4940 0.1840 0.5560]);
hold on
plot([vfBasket(fskB),vfBasket(fskB)],[fsC66Basket,fsC66Basket]./E,'Marker','.','MarkerSize',8,'LineStyle','none','Color',[0.4660 0.6740 0.1880]);
plot([vfBasket(askB),vfBasket(askB)],[C66Basket,C66Basket]./E,'Marker','.','MarkerSize',8,'LineStyle','none','Color',[0.3010 0.7450 0.9330]);
grid on
xlabel('Volume Fraction');
ylabel('Shear Stiffness');
legend('Double-Mirror Symmetric Designs','Mirror Symmetric Designs','Asymmetric Designs','Location','northeast');
axis([0.1 0.6 0 0.1]);
set(gca,'fontsize', 13)

%% PLOT RESULTS- 2D SCATTERPLOT (Poissons-vF)
plot([vfBasket(skB),vfBasket(skB)],[svyxBasket,svxyBasket],'Marker','.','MarkerSize',8,'LineStyle','none','Color',[0.4940 0.1840 0.5560]);
hold on
plot([vfBasket(fskB),vfBasket(fskB)],[fsvyxBasket,fsvxyBasket],'Marker','.','MarkerSize',8,'LineStyle','none','Color',[0.4660 0.6740 0.1880]);
plot([vfBasket(askB),vfBasket(askB)],[vyxBasket,vxyBasket],'Marker','.','MarkerSize',8,'LineStyle','none','Color',[0.3010 0.7450 0.9330]);
grid on
xlabel('Volume Fraction');
ylabel('Poisson''s Ratio');
legend('Double-Mirror Symmetric','Mirror Symmetric','Asymmetric','Location','southeast');
axis([0.1 0.6 -1 1]);
set(gca,'fontsize', 13)

%% PLOT SINGLE RESULTS- 2D SCATTERPLOTS (Poissons-Poissons)
plot([vxyBasket,vyxBasket],[vyxBasket,vxyBasket],'Marker','.','LineStyle','none','MarkerSize',8,'Color',[0.3010 0.7450 0.9330]);
hold on
plot([fsvxyBasket,fsvyxBasket],[fsvyxBasket,fsvxyBasket],'Marker','.','LineStyle','none','MarkerSize',8,'Color',[0.4660 0.6740 0.1880]);
plot([svxyBasket,svyxBasket],[svyxBasket,svxyBasket],'Marker','.','LineStyle','none','MarkerSize',8,'Color',[0.4940 0.1840 0.5560]);
grid on
xlabel('\nu_{12}');
ylabel('\nu_{21}');
legend('Asymmetric','Mirror Symmetric','Double-Mirror Symmetric','Location','northwest');
set(gca,'fontsize', 13)

%% PLOT SINGLE RESULTS- 2D SCATTERPLOTS (Normal-Normal)
plot([C11Basket,C22Basket]./E,[C22Basket,C11Basket]./E,'Marker','.','MarkerSize',8,'LineStyle','none','Color',[0.3010 0.7450 0.9330]);
hold on
plot([fsC11Basket,fsC22Basket]./E,[fsC22Basket,fsC11Basket]./E,'Marker','.','MarkerSize',8,'LineStyle','none','Color',[0.4660 0.6740 0.1880]);
plot([sC11Basket,sC22Basket]./E,[sC22Basket,sC11Basket]./E,'Marker','.','MarkerSize',8,'LineStyle','none','Color',[0.4940 0.1840 0.5560]);
grid on
xlabel('Horizontal Normal Stiffness');
ylabel('Vertical Normal Stiffness');
legend('Asymmetric','Mirror Symmetric','Double-Mirror Symmetric','Location','northwest');
set(gca,'fontsize', 13)
axis([0 0.4 0 0.4]);

%% PLOT SINGLE RESULTS- 2D SCATTERPLOTS (Normal-VF)
CaxBasket = [C11Basket,C22Basket]./E; 
sCaxBasket = [sC11Basket,sC22Basket]./E;
fsCaxBasket = [fsC11Basket,fsC22Basket]./E;

% Develop regressions
asfit = polyfit([vfBasket(askB),vfBasket(askB)],CaxBasket,1);          
asyfit = polyval(asfit,[vfBasket(askB),vfBasket(askB)]);          
asSStot = sum((CaxBasket-mean(CaxBasket)).^2);                   
asSSres = sum((CaxBasket-asyfit).^2);                       
asRsq = 1-asSSres/asSStot;          

fsfit = polyfit([vfBasket(fskB),vfBasket(fskB)],fsCaxBasket,1);          
fsyfit = polyval(fsfit,[vfBasket(fskB),vfBasket(fskB)]);          
fsSStot = sum((fsCaxBasket-mean(fsCaxBasket)).^2);                   
fsSSres = sum((fsCaxBasket-fsyfit).^2);                       
fsRsq = 1-fsSSres/fsSStot;   

sfit = polyfit([vfBasket(skB),vfBasket(skB)],sCaxBasket,1);          
syfit = polyval(sfit,[vfBasket(skB),vfBasket(skB)]);          
sSStot = sum((sCaxBasket-mean(sCaxBasket)).^2);                   
sSSres = sum((sCaxBasket-syfit).^2);                       
sRsq = 1-sSSres/sSStot; 

% Plotting
plot([vfBasket(askB),vfBasket(askB)],CaxBasket,'Marker','.','MarkerSize',8,'LineStyle','none','Color',[0.3010 0.7450 0.9330]);
hold on
plot([vfBasket(fskB),vfBasket(fskB)],fsCaxBasket,'Marker','.','MarkerSize',8,'LineStyle','none','Color',[0.4660 0.6740 0.1880]);
plot([vfBasket(skB),vfBasket(skB)],sCaxBasket,'Marker','.','MarkerSize',8,'LineStyle','none','Color',[0.4940 0.1840 0.5560]);
grid on
ylabel('Normal Stiffness');
xlabel('Volume Fraction');
legend('Asymmetric','Mirror Symmetric','Double-Mirror Symmetric','Location','northwest');
set(gca,'fontsize', 13)
axis([0.1 0.6 0 0.4]);
hold off
fprintf('[asymm slope, asymm R^2] = [%.5f, %.5f]\n',[asfit(1),asRsq])
fprintf('[fsymm slope, fsymm R^2] = [%.5f, %.5f]\n',[fsfit(1),fsRsq])
fprintf('[symm slope, symm R^2] = [%.5f, %.5f]\n',[sfit(1),sRsq])

%% PLOT SINGLE RESULTS- 2D SCATTERPLOTS (Normal-Poissons)
% Initialize plotting vectors
CaxBasket = [C11Basket,C22Basket]./E; 
sCaxBasket = [sC11Basket,sC22Basket]./E;
fsCaxBasket = [fsC11Basket,fsC22Basket]./E;
vbothBasket = [vyxBasket,vxyBasket];
svbothBasket = [svyxBasket,svxyBasket];
fsvbothBasket = [fsvyxBasket,fsvxyBasket];

% Plotting
plot(CaxBasket,vbothBasket,'Marker','.','MarkerSize',8,'LineStyle','none','Color',[0.3010 0.7450 0.9330]);
hold on
plot(fsCaxBasket,fsvbothBasket,'Marker','.','MarkerSize',8,'LineStyle','none','Color',[0.4660 0.6740 0.1880]);
plot(sCaxBasket,svbothBasket,'Marker','.','MarkerSize',8,'LineStyle','none','Color',[0.4940 0.1840 0.5560]);
grid on
xlabel('Normal Stiffness');
ylabel('Poisson''s Ratio');
legend('Asymmetric','Mirror Symmetric','Double-Mirror Symmetric');
set(gca,'fontsize', 13)
axis([0 0.4 -1 1]);
fprintf('[asymm mean, asymm min, asymm max] = [%.5f, %.5f, %.5f]\n',[mean(vbothBasket),min(vbothBasket),max(vbothBasket)])
fprintf('[fsymm mean, fsymm min, fsymm max] = [%.5f, %.5f, %.5f]\n',[mean(fsvbothBasket),min(fsvbothBasket),max(fsvbothBasket)])
fprintf('[symm mean, symm min, symm max] = [%.5f, %.5f, %.5f]\n',[mean(svbothBasket),min(svbothBasket),max(svbothBasket)])

%% PLOT SINGLE RESULTS- 2D SCATTERPLOTS (Normal-Shear)
% Initialize plotting vectors
CaxBasket = [C11Basket,C22Basket]./E; 
sCaxBasket = [sC11Basket,sC22Basket]./E;
fsCaxBasket = [fsC11Basket,fsC22Basket]./E;
CshBasket = [C66Basket,C66Basket]./E;
sCshBasket = [sC66Basket,sC66Basket]./E;
fsCshBasket = [fsC66Basket,fsC66Basket]./E;

% Develop regressions
asfit = polyfit(CaxBasket,CshBasket,1);          
asyfit = polyval(asfit,CaxBasket);          
asSStot = sum((CshBasket-mean(CshBasket)).^2);                   
asSSres = sum((CshBasket-asyfit).^2);                       
asRsq = 1-asSSres/asSStot;          

fsfit = polyfit(fsCaxBasket,fsCshBasket,1);          
fsyfit = polyval(fsfit,fsCaxBasket);          
fsSStot = sum((fsCshBasket-mean(fsCshBasket)).^2);                   
fsSSres = sum((fsCshBasket-fsyfit).^2);                       
fsRsq = 1-fsSSres/fsSStot;   

sfit = polyfit(sCaxBasket,sCshBasket,1);          
syfit = polyval(sfit,sCaxBasket);          
sSStot = sum((sCshBasket-mean(sCshBasket)).^2);                   
sSSres = sum((sCshBasket-syfit).^2);                       
sRsq = 1-sSSres/sSStot;   

% Plotting
plot(CaxBasket,CshBasket,'Marker','.','MarkerSize',8,'LineStyle','none','Color',[0.3010 0.7450 0.9330]);
hold on
plot(fsCaxBasket,fsCshBasket,'Marker','.','MarkerSize',8,'LineStyle','none','Color',[0.4660 0.6740 0.1880]);
plot(sCaxBasket,sCshBasket,'Marker','.','MarkerSize',8,'LineStyle','none','Color',[0.4940 0.1840 0.5560]);
grid on
xlabel('Normal Stiffness');
ylabel('Shear Stiffness');
legend('Asymmetric','Mirror Symmetric','Double-Mirror Symmetric');
set(gca,'fontsize', 13)
axis([0 0.4 0 0.08]);
fprintf('[asymm slope, asymm R^2] = [%.5f, %.5f]\n',[asfit(1),asRsq])
fprintf('[fsymm slope, fsymm R^2] = [%.5f, %.5f]\n',[fsfit(1),fsRsq])
fprintf('[symm slope, symm R^2] = [%.5f, %.5f]\n',[sfit(1),sRsq])

%% PDF Plots for Aggregated Values
fig = figure(1);

% First plot: Shear stiffness values
asCvalBasket = (2.*C66Basket)./(C11Basket+C22Basket);
asCvalBasket = asCvalBasket((asCvalBasket > 0) & (asCvalBasket < 1));
fsCvalBasket = (2.*fsC66Basket)./(fsC11Basket+fsC22Basket);
fsCvalBasket = fsCvalBasket((fsCvalBasket > 0) & (fsCvalBasket < 1));
sCvalBasket = (2.*sC66Basket)./(sC11Basket+sC22Basket);
sCvalBasket = sCvalBasket((sCvalBasket > 0) & (sCvalBasket < 1));

subplot(3,1,1)
x = 0:0.001:1;
[~,~,h_as] = ksdensity(([C11Basket,C22Basket]./E));
pdas = fitdist(([C11Basket,C22Basket]./E)','Kernel','Width',h_as);
[~,~,h_fs] = ksdensity(([fsC11Basket,fsC22Basket]./E));
pdfs = fitdist(([fsC11Basket,fsC22Basket]./E)','Kernel','Width',h_fs);
[~,~,h_s] = ksdensity(([sC11Basket,sC22Basket]./E));
pds = fitdist(([sC11Basket,sC22Basket]./E)','Kernel','Width',h_s);
yas = pdf(pdas,x);
yfs = pdf(pdfs,x);
ys = pdf(pds,x);
plot(x,yas,'LineStyle','-','Color',[0.3010 0.7450 0.9330],'LineWidth',2)
hold on
plot(x,yfs,'LineStyle','-','Color',[0.4660 0.6740 0.1880],'LineWidth',2)
plot(x,ys,'LineStyle','-','Color',[0.4940 0.1840 0.5560],'LineWidth',2)
xlabel('Normal Stiffness')
xlim([0 0.4]);
set(gca,'fontsize', 13)

% Second plot: Poisson's ratio
subplot(3,1,2)
x = -1:0.001:1;
[~,~,h_as] = ksdensity(vyxBasket);
pdas = fitdist(vyxBasket','Kernel','Width',h_as);
[~,~,h_fs] = ksdensity(fsvyxBasket);
pdfs = fitdist(fsvyxBasket','Kernel','Width',h_fs);
[~,~,h_s] = ksdensity(svyxBasket);
pds = fitdist(svyxBasket','Kernel','Width',h_s);
yas = pdf(pdas,x);
yfs = pdf(pdfs,x);
ys = pdf(pds,x);
plot(x,yas,'LineStyle','-','Color',[0.3010 0.7450 0.9330],'LineWidth',2)
hold on
plot(x,yfs,'LineStyle','-','Color',[0.4660 0.6740 0.1880],'LineWidth',2)
plot(x,ys,'LineStyle','-','Color',[0.4940 0.1840 0.5560],'LineWidth',2)
xlabel('Poisson''s Ratio')
xlim([-1 1]);
set(gca,'fontsize', 13)

% Third Plot: Normal Stiffness
subplot(3,1,3)
x = 0:0.001:0.3;
[~,~,h_as] = ksdensity((asCvalBasket));
pdas = fitdist((asCvalBasket)','Kernel','Width',h_as);
[~,~,h_fs] = ksdensity((fsCvalBasket));
pdfs = fitdist((fsCvalBasket)','Kernel','Width',h_fs);
[~,~,h_s] = ksdensity((sCvalBasket));
pds = fitdist((asCvalBasket)','Kernel','Width',h_s);
yas = pdf(pdas,x);
yfs = pdf(pdfs,x);
ys = pdf(pds,x);
plot(x,yas,'LineStyle','-','Color',[0.3010 0.7450 0.9330],'LineWidth',2)
hold on
plot(x,yfs,'LineStyle','-','Color',[0.4660 0.6740 0.1880],'LineWidth',2)
plot(x,ys,'LineStyle','-','Color',[0.4940 0.1840 0.5560],'LineWidth',2)
xlabel('Normalized Shear Stiffness')
xlim([0 0.3]);
set(gca,'fontsize', 13)

% Plot y-label
han=axes(fig,'visible','off'); 
han.YLabel.Visible='on';
ylabel(han,'Probability Density');
set(gca,'fontsize', 13)

% K-S tests 
aggregateKSRes = zeros(3,3);
aggregateKSRes(1,1) = kstest2(fsCvalBasket,asCvalBasket);
aggregateKSRes(1,2) = kstest2(fsCvalBasket,sCvalBasket);
aggregateKSRes(1,3) = kstest2(asCvalBasket,sCvalBasket);
aggregateKSRes(2,1) = kstest2(vyxBasket,fsvyxBasket);
aggregateKSRes(2,2) = kstest2(svyxBasket,vyxBasket);
aggregateKSRes(2,3) = kstest2(svyxBasket,fsvyxBasket);
aggregateKSRes(3,1) = kstest2([C11Basket,C22Basket]./E,[fsC11Basket,fsC22Basket]./E);
aggregateKSRes(3,2) = kstest2([C11Basket,C22Basket]./E,[sC11Basket,sC22Basket]./E);
aggregateKSRes(3,3) = kstest2([fsC11Basket,fsC22Basket]./E,[sC11Basket,sC22Basket]./E);

%% Design Characteristics Finder (Binary)
% This script examines all fully-symmetric, flip-symmetric, and
% asymmetric designs within a design space for the number of key
% characteristics associated with symmetry

% Initialize values
desiredChars = [2,4,5,6];

allC66Basket = [C66Basket,sC66Basket,fsC66Basket];
allC26Basket = [C26Basket,sC26Basket,fsC26Basket];
allC22Basket = [C22Basket,sC22Basket,fsC22Basket];
allC11Basket = [C11Basket,sC11Basket,fsC11Basket];
allvyxBasket = [vyxBasket,svyxBasket,fsvyxBasket];
allCABasket = [goodCABasket,sgoodCABasket,fsgoodCABasket];
allCABasket = allCABasket(~cellfun('isempty',allCABasket));
usedCAs = {}; uCAc = 1;

asCharCountBins = []; fsCharCountBins = []; sCharCountBins = [];
uasc = 0; usc = 0; ufsc = 0;
uniqueASCA = {}; uniqueFSCA = {}; uniqueSCA = {};
uniqueASidx = []; uniqueFSidx = []; uniqueSidx = [];

% Loop through all designs
for v = 1:1:length(allCABasket)
    tic
    % Initialize values
    CA = cell2mat(allCABasket(v));
    
    % Check if current CA has been used already
    usedbool = 0;
    for k = 1:1:length(usedCAs)
        uCA = cell2mat(usedCAs(k));
        if size(CA,1) == size(uCA,1)
            if CA == uCA
                usedbool = 1;
                break
            end
        end
    end
    if usedbool == 0
        % Evaluate booleans
        charCounts = desCharFinder_3x3(CA,NC,sel,sidenum,desiredChars);
        % Evaluate symmetry score
        symmScore = symmHeuristic_2D(CA,NC,sel);
        % Record boolean results
        if symmScore == 0
            asCharCountBins = [asCharCountBins;charCounts];
            uasc = uasc + 1;
            uniqueASCA(uasc) = {CA};
            uniqueASidx = [uniqueASidx,v];
        elseif symmScore == 1
            fsCharCountBins = [fsCharCountBins;charCounts];
            ufsc = ufsc + 1;
            uniqueFSCA(ufsc) = {CA};
            uniqueFSidx = [uniqueFSidx,v];
        elseif symmScore == 2
            sCharCountBins = [sCharCountBins;charCounts];
            usc = usc + 1;
            uniqueSCA(usc) = {CA};
            uniqueSidx = [uniqueSidx,v];
        end
    end
    
    % Log current CA as used
    usedCAs(uCAc) = {CA}; uCAc = uCAc + 1;
    disp('Current Design'); disp(v);
    toc
end

asymmCharCounts = [];
for i = 1:4:(uasc*4)
    row = asCharCountBins(i:(i+3));
    asymmCharCounts = [asymmCharCounts;row'];
end
fsymmCharCounts = [];
for i = 1:4:(ufsc*4)
    row = fsCharCountBins(i:(i+3));
    fsymmCharCounts = [fsymmCharCounts;row'];
end
symmCharCounts = [];
for i = 1:4:(usc*4)
    row = sCharCountBins(i:(i+3));
    symmCharCounts = [symmCharCounts;row'];
end
                
%% Plotting Design Characteristics Bar Chart

% Randomize indices of each data set
numParts = 5;
asymmScrambleIdx = randperm(uasc); 
asSplitSize = round((uasc./numParts),0);
fsymmScrambleIdx = randperm(ufsc); 
fsSplitSize = round((ufsc./numParts),0);
symmScrambleIdx = randperm(usc); 
sSplitSize = round((usc./numParts),0);
asNormVals = []; fsNormVals = []; sNormVals = [];

% Find the normalized proportion of designs with each characteristic in
% each subset
for i = 1:1:(numParts-1)
    asymmVals = asymmCharCounts((((i-1)*asSplitSize)+1):(i*asSplitSize),:);
    asRow = mean(asymmVals);
    fsymmVals = fsymmCharCounts((((i-1)*fsSplitSize)+1):(i*fsSplitSize),:);
    fsRow = mean(fsymmVals);
    symmVals = symmCharCounts((((i-1)*sSplitSize)+1):(i*sSplitSize),:);
    sRow = mean(symmVals);
    asNormVals = [asNormVals;asRow];
    fsNormVals = [fsNormVals;fsRow];
    sNormVals = [sNormVals;sRow];
end
finalasVals = asymmCharCounts((((numParts-1)*asSplitSize)+1):end,:);
finalasRow = mean(finalasVals);
finalfsVals = fsymmCharCounts((((numParts-1)*fsSplitSize)+1):end,:);
finalfsRow = mean(finalfsVals);
finalsVals = symmCharCounts((((numParts-1)*sSplitSize)+1):end,:);
finalsRow = mean(finalsVals);
asNormVals = [asNormVals;finalasRow];
fsNormVals = [fsNormVals;finalfsRow];
sNormVals = [sNormVals;finalsRow];

% Find the means and 95% error of each set
asMeans = mean(asNormVals); as95p = sqrt(var(asNormVals));
fsMeans = mean(fsNormVals); fs95p = sqrt(var(fsNormVals));
sMeans = mean(sNormVals); s95p = sqrt(var(sNormVals));
data = [asMeans',fsMeans',sMeans'];
error = [as95p;fs95p;s95p];

% Plot design characteristics frequency
barlabs = categorical({'Pivot Points','Arrow Shapes','Spider Nodes','Stacked Members'});
barlabs = reordercats(barlabs,{'Pivot Points','Arrow Shapes','Spider Nodes','Stacked Members'});
b = bar(barlabs,data);
b(1).FaceColor = [0.3010 0.7450 0.9330];
b(2).FaceColor = [0.4660 0.6740 0.1880];
b(3).FaceColor = [0.4940 0.1840 0.5560];
hBar = bar(data, 0.8);                                                    
for k1 = 1:size(data,2)
    ctr(k1,:) = bsxfun(@plus, hBar(k1).XData, hBar(k1).XOffset');       
    ydt(k1,:) = hBar(k1).YData;                                       
end
hold on
errorbar(ctr,ydt,error,'.k');
legend('Asymmetry','Mirror Symmetry','Double-Mirror Symmetry','Location','northwest');
ylabel('Normalized Frequency');
set(gca,'xticklabel',barlabs);
hBar(1).FaceColor = [0.3010 0.7450 0.9330];
hBar(2).FaceColor = [0.4660 0.6740 0.1880];
hBar(3).FaceColor = [0.4940 0.1840 0.5560];
set(gca,'fontsize', 13)
hold off

%% Formatting Binary Design Characteristics Data

% Characteristic 2
asC66_with_2 = asymmCharCounts(:,1).*allC66Basket(uniqueASidx)';
asC22_with_2 = (asymmCharCounts(:,1).*allC22Basket(uniqueASidx)');
asC11_with_2 = (asymmCharCounts(:,1).*allC11Basket(uniqueASidx)');
asvyx_with_2 = (asymmCharCounts(:,1).*allvyxBasket(uniqueASidx)');
fsC66_with_2 = (fsymmCharCounts(:,1).*allC66Basket(uniqueFSidx)');
fsC22_with_2 = (fsymmCharCounts(:,1).*allC22Basket(uniqueFSidx)');
fsC11_with_2 = (fsymmCharCounts(:,1).*allC11Basket(uniqueFSidx)');
fsvyx_with_2 = (fsymmCharCounts(:,1).*allvyxBasket(uniqueFSidx)');
sC66_with_2 = (symmCharCounts(:,1).*allC66Basket(uniqueSidx)');
sC22_with_2 = (symmCharCounts(:,1).*allC22Basket(uniqueSidx)');
sC11_with_2 = (symmCharCounts(:,1).*allC11Basket(uniqueSidx)');
svyx_with_2 = (symmCharCounts(:,1).*allvyxBasket(uniqueSidx)');

asC66_with_2 = asC66_with_2(asC66_with_2 ~= 0);
asC22_with_2 = asC22_with_2(asC22_with_2 ~= 0);
asC11_with_2 = asC11_with_2(asC11_with_2 ~= 0);
asvyx_with_2 = asvyx_with_2(asvyx_with_2 ~= 0);
fsC66_with_2 = fsC66_with_2(fsC66_with_2 ~= 0);
fsC22_with_2 = fsC22_with_2(fsC22_with_2 ~= 0);
fsC11_with_2 = fsC11_with_2(fsC11_with_2 ~= 0);
fsvyx_with_2 = fsvyx_with_2(fsvyx_with_2 ~= 0);
sC66_with_2 = sC66_with_2(sC66_with_2 ~= 0);
sC22_with_2 = sC22_with_2(sC22_with_2 ~= 0);
sC11_with_2 = sC11_with_2(sC11_with_2 ~= 0);
svyx_with_2 = svyx_with_2(svyx_with_2 ~= 0);

asC66_without_2 = ((-1).*(asymmCharCounts(:,1)-1).*allC66Basket(uniqueASidx)');
asC22_without_2 = ((-1).*(asymmCharCounts(:,1)-1).*allC22Basket(uniqueASidx)');
asC11_without_2 = ((-1).*(asymmCharCounts(:,1)-1).*allC11Basket(uniqueASidx)');
asvyx_without_2 = ((-1).*(asymmCharCounts(:,1)-1).*allvyxBasket(uniqueASidx)');
fsC66_without_2 = ((-1).*(fsymmCharCounts(:,1)-1).*allC66Basket(uniqueFSidx)');
fsC22_without_2 = ((-1).*(fsymmCharCounts(:,1)-1).*allC22Basket(uniqueFSidx)');
fsC11_without_2 = ((-1).*(fsymmCharCounts(:,1)-1).*allC11Basket(uniqueFSidx)');
fsvyx_without_2 = ((-1).*(fsymmCharCounts(:,1)-1).*allvyxBasket(uniqueFSidx)');
sC66_without_2 = ((-1).*(symmCharCounts(:,1)-1).*allC66Basket(uniqueSidx)');
sC22_without_2 = ((-1).*(symmCharCounts(:,1)-1).*allC22Basket(uniqueSidx)');
sC11_without_2 = ((-1).*(symmCharCounts(:,1)-1).*allC11Basket(uniqueSidx)');
svyx_without_2 = ((-1).*(symmCharCounts(:,1)-1).*allvyxBasket(uniqueSidx)');

asC66_without_2 = asC66_without_2(asC66_without_2 ~= 0);
asC22_without_2 = asC22_without_2(asC22_without_2 ~= 0);
asC11_without_2 = asC11_without_2(asC11_without_2 ~= 0);
asvyx_without_2 = asvyx_without_2(asvyx_without_2 ~= 0);
fsC66_without_2 = fsC66_without_2(fsC66_without_2 ~= 0);
fsC22_without_2 = fsC22_without_2(fsC22_without_2 ~= 0);
fsC11_without_2 = fsC11_without_2(fsC11_without_2 ~= 0);
fsvyx_without_2 = fsvyx_without_2(fsvyx_without_2 ~= 0);
sC66_without_2 = sC66_without_2(sC66_without_2 ~= 0);
sC22_without_2 = sC22_without_2(sC22_without_2 ~= 0);
sC11_without_2 = sC11_without_2(sC11_without_2 ~= 0);
svyx_without_2 = svyx_without_2(svyx_without_2 ~= 0);

% Characteristic 4
asC66_with_4 = (asymmCharCounts(:,2).*allC66Basket(uniqueASidx)');
asC22_with_4 = (asymmCharCounts(:,2).*allC22Basket(uniqueASidx)');
asC11_with_4 = (asymmCharCounts(:,2).*allC11Basket(uniqueASidx)');
asvyx_with_4 = (asymmCharCounts(:,2).*allvyxBasket(uniqueASidx)');
fsC66_with_4 = (fsymmCharCounts(:,2).*allC66Basket(uniqueFSidx)');
fsC22_with_4 = (fsymmCharCounts(:,2).*allC22Basket(uniqueFSidx)');
fsC11_with_4 = (fsymmCharCounts(:,2).*allC11Basket(uniqueFSidx)');
fsvyx_with_4 = (fsymmCharCounts(:,2).*allvyxBasket(uniqueFSidx)');
sC66_with_4 = (symmCharCounts(:,2).*allC66Basket(uniqueSidx)');
sC22_with_4 = (symmCharCounts(:,2).*allC22Basket(uniqueSidx)');
sC11_with_4 = (symmCharCounts(:,2).*allC11Basket(uniqueSidx)');
svyx_with_4 = (symmCharCounts(:,2).*allvyxBasket(uniqueSidx)');

asC66_with_4 = asC66_with_4(asC66_with_4 ~= 0);
asC22_with_4 = asC22_with_4(asC22_with_4 ~= 0);
asC11_with_4 = asC11_with_4(asC11_with_4 ~= 0);
asvyx_with_4 = asvyx_with_4(asvyx_with_4 ~= 0);
fsC66_with_4 = fsC66_with_4(fsC66_with_4 ~= 0);
fsC22_with_4 = fsC22_with_4(fsC22_with_4 ~= 0);
fsC11_with_4 = fsC11_with_4(fsC11_with_4 ~= 0);
fsvyx_with_4 = fsvyx_with_4(fsvyx_with_4 ~= 0);
sC66_with_4 = sC66_with_4(sC66_with_4 ~= 0);
sC22_with_4 = sC22_with_4(sC22_with_4 ~= 0);
sC11_with_4 = sC11_with_4(sC11_with_4 ~= 0);
svyx_with_4 = svyx_with_4(svyx_with_4 ~= 0);

asC66_without_4 = ((-1).*(asymmCharCounts(:,2)-1).*allC66Basket(uniqueASidx)');
asC22_without_4 = ((-1).*(asymmCharCounts(:,2)-1).*allC22Basket(uniqueASidx)');
asC11_without_4 = ((-1).*(asymmCharCounts(:,2)-1).*allC11Basket(uniqueASidx)');
asvyx_without_4 = ((-1).*(asymmCharCounts(:,2)-1).*allvyxBasket(uniqueASidx)');
fsC66_without_4 = ((-1).*(fsymmCharCounts(:,2)-1).*allC66Basket(uniqueFSidx)');
fsC22_without_4 = ((-1).*(fsymmCharCounts(:,2)-1).*allC22Basket(uniqueFSidx)');
fsC11_without_4 = ((-1).*(fsymmCharCounts(:,2)-1).*allC11Basket(uniqueFSidx)');
fsvyx_without_4 = ((-1).*(fsymmCharCounts(:,2)-1).*allvyxBasket(uniqueFSidx)');
sC66_without_4 = ((-1).*(symmCharCounts(:,2)-1).*allC66Basket(uniqueSidx)');
sC22_without_4 = ((-1).*(symmCharCounts(:,2)-1).*allC22Basket(uniqueSidx)');
sC11_without_4 = ((-1).*(symmCharCounts(:,2)-1).*allC11Basket(uniqueSidx)');
svyx_without_4 = ((-1).*(symmCharCounts(:,2)-1).*allvyxBasket(uniqueSidx)');

asC66_without_4 = asC66_without_4(asC66_without_4 ~= 0);
asC22_without_4 = asC22_without_4(asC22_without_4 ~= 0);
asC11_without_4 = asC11_without_4(asC11_without_4 ~= 0);
asvyx_without_4 = asvyx_without_4(asvyx_without_4 ~= 0);
fsC66_without_4 = fsC66_without_4(fsC66_without_4 ~= 0);
fsC22_without_4 = fsC22_without_4(fsC22_without_4 ~= 0);
fsC11_without_4 = fsC11_without_4(fsC11_without_4 ~= 0);
fsvyx_without_4 = fsvyx_without_4(fsvyx_without_4 ~= 0);
sC66_without_4 = sC66_without_4(sC66_without_4 ~= 0);
sC22_without_4 = sC22_without_4(sC22_without_4 ~= 0);
sC11_without_4 = sC11_without_4(sC11_without_4 ~= 0);
svyx_without_4 = svyx_without_4(svyx_without_4 ~= 0);

% Characteristic 5
asC66_with_5 = (asymmCharCounts(:,3).*allC66Basket(uniqueASidx)');
asC22_with_5 = (asymmCharCounts(:,3).*allC22Basket(uniqueASidx)');
asC11_with_5 = (asymmCharCounts(:,3).*allC11Basket(uniqueASidx)');
asvyx_with_5 = (asymmCharCounts(:,3).*allvyxBasket(uniqueASidx)');
fsC66_with_5 = (fsymmCharCounts(:,3).*allC66Basket(uniqueFSidx)');
fsC22_with_5 = (fsymmCharCounts(:,3).*allC22Basket(uniqueFSidx)');
fsC11_with_5 = (fsymmCharCounts(:,3).*allC11Basket(uniqueFSidx)');
fsvyx_with_5 = (fsymmCharCounts(:,3).*allvyxBasket(uniqueFSidx)');
sC66_with_5 = (symmCharCounts(:,3).*allC66Basket(uniqueSidx)');
sC22_with_5 = (symmCharCounts(:,3).*allC22Basket(uniqueSidx)');
sC11_with_5 = (symmCharCounts(:,3).*allC11Basket(uniqueSidx)');
svyx_with_5 = (symmCharCounts(:,3).*allvyxBasket(uniqueSidx)');

asC66_with_5 = asC66_with_5(asC66_with_5 ~= 0);
asC22_with_5 = asC22_with_5(asC22_with_5 ~= 0);
asC11_with_5 = asC11_with_5(asC11_with_5 ~= 0);
asvyx_with_5 = asvyx_with_5(asvyx_with_5 ~= 0);
fsC66_with_5 = fsC66_with_5(fsC66_with_5 ~= 0);
fsC22_with_5 = fsC22_with_5(fsC22_with_5 ~= 0);
fsC11_with_5 = fsC11_with_5(fsC11_with_5 ~= 0);
fsvyx_with_5 = fsvyx_with_5(fsvyx_with_5 ~= 0);
sC66_with_5 = sC66_with_5(sC66_with_5 ~= 0);
sC22_with_5 = sC22_with_5(sC22_with_5 ~= 0);
sC11_with_5 = sC11_with_5(sC11_with_5 ~= 0);
svyx_with_5 = svyx_with_5(svyx_with_5 ~= 0);

asC66_without_5 = ((-1).*(asymmCharCounts(:,3)-1).*allC66Basket(uniqueASidx)');
asC22_without_5 = ((-1).*(asymmCharCounts(:,3)-1).*allC22Basket(uniqueASidx)');
asC11_without_5 = ((-1).*(asymmCharCounts(:,3)-1).*allC11Basket(uniqueASidx)');
asvyx_without_5 = ((-1).*(asymmCharCounts(:,3)-1).*allvyxBasket(uniqueASidx)');
fsC66_without_5 = ((-1).*(fsymmCharCounts(:,3)-1).*allC66Basket(uniqueFSidx)');
fsC22_without_5 = ((-1).*(fsymmCharCounts(:,3)-1).*allC22Basket(uniqueFSidx)');
fsC11_without_5 = ((-1).*(fsymmCharCounts(:,3)-1).*allC11Basket(uniqueFSidx)');
fsvyx_without_5 = ((-1).*(fsymmCharCounts(:,3)-1).*allvyxBasket(uniqueFSidx)');
sC66_without_5 = ((-1).*(symmCharCounts(:,3)-1).*allC66Basket(uniqueSidx)');
sC22_without_5 = ((-1).*(symmCharCounts(:,3)-1).*allC22Basket(uniqueSidx)');
sC11_without_5 = ((-1).*(symmCharCounts(:,3)-1).*allC11Basket(uniqueSidx)');
svyx_without_5 = ((-1).*(symmCharCounts(:,3)-1).*allvyxBasket(uniqueSidx)');

asC66_without_5 = asC66_without_5(asC66_without_5 ~= 0);
asC22_without_5 = asC22_without_5(asC22_without_5 ~= 0);
asC11_without_5 = asC11_without_5(asC11_without_5 ~= 0);
asvyx_without_5 = asvyx_without_5(asvyx_without_5 ~= 0);
fsC66_without_5 = fsC66_without_5(fsC66_without_5 ~= 0);
fsC22_without_5 = fsC22_without_5(fsC22_without_5 ~= 0);
fsC11_without_5 = fsC11_without_5(fsC11_without_5 ~= 0);
fsvyx_without_5 = fsvyx_without_5(fsvyx_without_5 ~= 0);
sC66_without_5 = sC66_without_5(sC66_without_5 ~= 0);
sC22_without_5 = sC22_without_5(sC22_without_5 ~= 0);
sC11_without_5 = sC11_without_5(sC11_without_5 ~= 0);
svyx_without_5 = svyx_without_5(svyx_without_5 ~= 0);

% Characteristic 6
asC66_with_6 = (asymmCharCounts(:,4).*allC66Basket(uniqueASidx)');
asC22_with_6 = (asymmCharCounts(:,4).*allC22Basket(uniqueASidx)');
asC11_with_6 = (asymmCharCounts(:,4).*allC11Basket(uniqueASidx)');
asvyx_with_6 = (asymmCharCounts(:,4).*allvyxBasket(uniqueASidx)');
fsC66_with_6 = (fsymmCharCounts(:,4).*allC66Basket(uniqueFSidx)');
fsC22_with_6 = (fsymmCharCounts(:,4).*allC22Basket(uniqueFSidx)');
fsC11_with_6 = (fsymmCharCounts(:,4).*allC11Basket(uniqueFSidx)');
fsvyx_with_6 = (fsymmCharCounts(:,4).*allvyxBasket(uniqueFSidx)');
sC66_with_6 = (symmCharCounts(:,4).*allC66Basket(uniqueSidx)');
sC22_with_6 = (symmCharCounts(:,4).*allC22Basket(uniqueSidx)');
sC11_with_6 = (symmCharCounts(:,4).*allC11Basket(uniqueSidx)');
svyx_with_6 = (symmCharCounts(:,4).*allvyxBasket(uniqueSidx)');

asC66_with_6 = asC66_with_6(asC66_with_6 ~= 0);
asC22_with_6 = asC22_with_6(asC22_with_6 ~= 0);
asC11_with_6 = asC11_with_6(asC11_with_6 ~= 0);
asvyx_with_6 = asvyx_with_6(asvyx_with_6 ~= 0);
fsC66_with_6 = fsC66_with_6(fsC66_with_6 ~= 0);
fsC22_with_6 = fsC22_with_6(fsC22_with_6 ~= 0);
fsC11_with_6 = fsC11_with_6(fsC11_with_6 ~= 0);
fsvyx_with_6 = fsvyx_with_6(fsvyx_with_6 ~= 0);
sC66_with_6 = sC66_with_6(sC66_with_6 ~= 0);
sC22_with_6 = sC22_with_6(sC22_with_6 ~= 0);
sC11_with_6 = sC11_with_6(sC11_with_6 ~= 0);
svyx_with_6 = svyx_with_6(svyx_with_6 ~= 0);

asC66_without_6 = ((-1).*(asymmCharCounts(:,4)-1).*allC66Basket(uniqueASidx)');
asC22_without_6 = ((-1).*(asymmCharCounts(:,4)-1).*allC22Basket(uniqueASidx)');
asC11_without_6 = ((-1).*(asymmCharCounts(:,4)-1).*allC11Basket(uniqueASidx)');
asvyx_without_6 = ((-1).*(asymmCharCounts(:,4)-1).*allvyxBasket(uniqueASidx)');
fsC66_without_6 = ((-1).*(fsymmCharCounts(:,4)-1).*allC66Basket(uniqueFSidx)');
fsC22_without_6 = ((-1).*(fsymmCharCounts(:,4)-1).*allC22Basket(uniqueFSidx)');
fsC11_without_6 = ((-1).*(fsymmCharCounts(:,4)-1).*allC11Basket(uniqueFSidx)');
fsvyx_without_6 = ((-1).*(fsymmCharCounts(:,4)-1).*allvyxBasket(uniqueFSidx)');
sC66_without_6 = ((-1).*(symmCharCounts(:,4)-1).*allC66Basket(uniqueSidx)');
sC22_without_6 = ((-1).*(symmCharCounts(:,4)-1).*allC22Basket(uniqueSidx)');
sC11_without_6 = ((-1).*(symmCharCounts(:,4)-1).*allC11Basket(uniqueSidx)');
svyx_without_6 = ((-1).*(symmCharCounts(:,4)-1).*allvyxBasket(uniqueSidx)');

asC66_without_6 = asC66_without_6(asC66_without_6 ~= 0);
asC22_without_6 = asC22_without_6(asC22_without_6 ~= 0);
asC11_without_6 = asC11_without_6(asC11_without_6 ~= 0);
asvyx_without_6 = asvyx_without_6(asvyx_without_6 ~= 0);
fsC66_without_6 = fsC66_without_6(fsC66_without_6 ~= 0);
fsC22_without_6 = fsC22_without_6(fsC22_without_6 ~= 0);
fsC11_without_6 = fsC11_without_6(fsC11_without_6 ~= 0);
fsvyx_without_6 = fsvyx_without_6(fsvyx_without_6 ~= 0);
sC66_without_6 = sC66_without_6(sC66_without_6 ~= 0);
sC22_without_6 = sC22_without_6(sC22_without_6 ~= 0);
sC11_without_6 = sC11_without_6(sC11_without_6 ~= 0);
svyx_without_6 = svyx_without_6(svyx_without_6 ~= 0);

%% Kernel-Fit Plots of Design Characteristics vs Properties (3x3)

%%% FIRST PLOT: Normalized C66 vs Pivot Points
fig1 = figure(1);

asCval_with_2 = (2.*asC66_with_2)./(asC11_with_2+asC22_with_2);
asCval_without_2 = (2.*asC66_without_2)./(asC11_without_2+asC22_without_2);

fsCval_with_2 = (2.*fsC66_with_2)./(fsC11_with_2+fsC22_with_2);
fsCval_without_2 = (2.*fsC66_without_2)./(fsC11_without_2+fsC22_without_2);

sCval_with_2 = (2.*sC66_with_2)./(sC11_with_2+sC22_with_2);
sCval_without_2 = (2.*sC66_without_2)./(sC11_without_2+sC22_without_2);

% First plot
subplot(3,1,1); 
x = 0:0.001:0.3;
[~,~,h_with] = ksdensity((asCval_with_2));
pdwith = fitdist((asCval_with_2),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity((asCval_without_2));
pdwithout = fitdist((asCval_without_2),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normalized Shear Stiffness')
legend('With Pivot Points','Without Pivot Points');
title('Asymmetric Designs');
xlim([0 0.3]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Second plot
subplot(3,1,2); 
x = 0:0.001:0.3;
[~,~,h_with] = ksdensity((fsCval_with_2));
pdwith = fitdist((fsCval_with_2),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity((fsCval_without_2));
pdwithout = fitdist((fsCval_without_2),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normalized Shear Stiffness')
title('Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Third plot
subplot(3,1,3); 
x = 0:0.001:0.3;
[~,~,h_with] = ksdensity((sCval_with_2));
pdwith = fitdist((sCval_with_2),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity((sCval_without_2));
pdwithout = fitdist((sCval_without_2),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normalized Shear Stiffness')
title('Double-Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'fontsize', 12)

han=axes(fig1,'visible','off'); 
han.YLabel.Visible='on';
ylabel(han,'Probability Density');
set(gca,'fontsize', 12)

%%% SECOND PLOT: v21 vs Pivot Points
fig2 = figure(2);

% First plot
subplot(3,1,1); 
x = -1:0.001:1;
[~,~,h_with] = ksdensity(asvyx_with_2);
pdwith = fitdist(asvyx_with_2,'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(asvyx_without_2);
pdwithout = fitdist(asvyx_without_2,'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Poisson''s Ratio')
legend('With Pivot Points','Without Pivot Points','Location','northwest');
title('Asymmetric Designs');
xlim([-0.6 0.6]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Second plot
subplot(3,1,2); 
x = -1:0.001:1;
[~,~,h_with] = ksdensity(fsvyx_with_2);
pdwith = fitdist(fsvyx_with_2,'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(fsvyx_without_2);
pdwithout = fitdist(fsvyx_without_2,'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
title('Mirror Symmetric Designs');
xlabel('Poisson''s Ratio')
xlim([-0.6 0.6]);
set(gca,'fontsize', 12)

% Third plot
subplot(3,1,3); 
x = -1:0.001:1;
[~,~,h_with] = ksdensity(svyx_with_2);
pdwith = fitdist(svyx_with_2,'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(svyx_without_2);
pdwithout = fitdist(svyx_without_2,'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
title('Double-Mirror Symmetric Designs');
xlabel('Poisson''s Ratio')
xlim([-0.6 0.6]);
set(gca,'fontsize', 12)

han=axes(fig2,'visible','off'); 
han.YLabel.Visible='on';
ylabel(han,'Probability Density');
set(gca,'fontsize', 12)

%%% THIRD PLOT: Normalized C66 vs Arrows
fig3 = figure(3);

asCval_with_4 = (2.*asC66_with_4)./(asC11_with_4+asC22_with_4);
asCval_without_4 = (2.*asC66_without_4)./(asC11_without_4+asC22_without_4);

fsCval_with_4 = (2.*fsC66_with_4)./(fsC11_with_4+fsC22_with_4);
fsCval_without_4 = (2.*fsC66_without_4)./(fsC11_without_4+fsC22_without_4);

sCval_with_4 = (2.*sC66_with_4)./(sC11_with_4+sC22_with_4);
sCval_without_4 = (2.*sC66_without_4)./(sC11_without_4+sC22_without_4);

% First plot
subplot(3,1,1); 
x = 0:0.001:0.3;
[~,~,h_with] = ksdensity((asCval_with_4));
pdwith = fitdist((asCval_with_4),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity((asCval_without_4));
pdwithout = fitdist((asCval_without_4),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
title('Asymmetric Designs');
xlabel('Normalized Shear Stiffness')
legend('With Arrows','Without Arrows');
xlim([0 0.3]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Second plot
subplot(3,1,2); 
x = 0:0.001:0.3;
[~,~,h_with] = ksdensity((fsCval_with_4));
pdwith = fitdist((fsCval_with_4),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity((fsCval_without_4));
pdwithout = fitdist((fsCval_without_4),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normalized Shear Stiffness')
title('Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Third plot
subplot(3,1,3); 
x = 0:0.001:0.3;
[~,~,h_without] = ksdensity((sCval_without_4));
pdwithout = fitdist((sCval_without_4),'Kernel','Width',h_without);
ywithout = pdf(pdwithout,x);
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normalized Shear Stiffness')
title('Double-Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'fontsize', 12)

han=axes(fig3,'visible','off'); 
han.YLabel.Visible='on';
ylabel(han,'Probability Density');
set(gca,'fontsize', 12)

%%% FOURTH PLOT: v21 vs Arrows
fig4 = figure(4);

% First plot
subplot(3,1,1); 
x = -1:0.001:1;
[~,~,h_with] = ksdensity(asvyx_with_4);
pdwith = fitdist(asvyx_with_4,'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(asvyx_without_4);
pdwithout = fitdist(asvyx_without_4,'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
title('Asymmetric Designs');
xlabel('Poisson''s Ratio');
legend('With Arrows','Without Arrows','Location','northwest');
xlim([-0.6 0.6]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Second plot
subplot(3,1,2); 
x = -1:0.001:1;
[~,~,h_with] = ksdensity(fsvyx_with_4);
pdwith = fitdist(fsvyx_with_4,'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(fsvyx_without_4);
pdwithout = fitdist(fsvyx_without_4,'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Poisson''s Ratio');
title('Mirror Symmetric Designs');
xlim([-0.6 0.6]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Third plot
subplot(3,1,3); 
x = -1:0.001:1;
[~,~,h_without] = ksdensity(svyx_without_4);
pdwithout = fitdist(svyx_without_4,'Kernel','Width',h_without);
ywithout = pdf(pdwithout,x);
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Poisson''s Ratio')
title('Double-Mirror Symmetric Designs');
xlim([-0.6 0.6]);
set(gca,'fontsize', 12)

han=axes(fig4,'visible','off'); 
han.YLabel.Visible='on';
ylabel(han,'Probability Density');
set(gca,'fontsize', 12)

%%% FIFTH PLOT: Normalized C66 vs Spider Nodes
fig5 = figure(5);

asCval_with_5 = (2.*asC66_with_5)./(asC11_with_5+asC22_with_5);
asCval_without_5 = (2.*asC66_without_5)./(asC11_without_5+asC22_without_5);

fsCval_with_5 = (2.*fsC66_with_5)./(fsC11_with_5+fsC22_with_5);
fsCval_without_5 = (2.*fsC66_without_5)./(fsC11_without_5+fsC22_without_5);

sCval_without_5 = (2.*sC66_without_5)./(sC11_without_5+sC22_without_5);

% First plot
subplot(3,1,1); 
x = 0:0.001:0.3;
[~,~,h_with] = ksdensity(rmoutliers(asCval_with_5));
pdwith = fitdist(rmoutliers(asCval_with_5),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(rmoutliers(asCval_without_5));
pdwithout = fitdist(rmoutliers(asCval_without_5),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
title('Asymmetric Designs');
xlabel('Normalized Shear Stiffness')
legend('With Spider Nodes','Without Spider Nodes');
xlim([0 0.3]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Second plot
subplot(3,1,2); 
x = 0:0.001:0.3;
[~,~,h_with] = ksdensity(rmoutliers(fsCval_with_5));
pdwith = fitdist(rmoutliers(fsCval_with_5),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(rmoutliers(fsCval_without_5));
pdwithout = fitdist(rmoutliers(fsCval_without_5),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normalized Shear Stiffness')
title('Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'fontsize', 12)

% Third plot
subplot(3,1,3); 
x = 0:0.001:0.3;
[~,~,h_without] = ksdensity((sCval_without_5));
pdwithout = fitdist((sCval_without_5),'Kernel','Width',h_without);
ywithout = pdf(pdwithout,x);
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normalized Shear Stiffness')
title('Double-Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'fontsize', 12)

han=axes(fig5,'visible','off'); 
han.YLabel.Visible='on';
ylabel(han,'Probability Density');
set(gca,'fontsize', 12)


%%% SIXTH PLOT: v21 vs Spider Nodes
fig6 = figure(6);

% First plot
subplot(3,1,1); 
x = -1:0.001:1;
pdwith = fitdist((asvyx_with_5),'Kernel','kernel','Normal');
pdwithout = fitdist((asvyx_without_5),'Kernel','kernel','Normal');
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
title('Asymmetric Designs');
xlabel('Poisson''s Ratio')
legend('With Spider Nodes','Without Spider Nodes','Location','northwest');
xlim([-0.6 0.6]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Second plot
subplot(3,1,2); 
x = -1:0.001:1;
[~,~,h_with] = ksdensity(fsvyx_with_5);
pdwith = fitdist(fsvyx_with_5,'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(fsvyx_without_5);
pdwithout = fitdist(fsvyx_without_5,'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Poisson''s Ratio')
title('Mirror Symmetric Designs');
xlim([-0.6 0.6]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Third plot
subplot(3,1,3); 
x = -1:0.001:1;
[~,~,h_without] = ksdensity(svyx_without_5);
pdwithout = fitdist(svyx_without_5,'Kernel','Width',h_without);
ywithout = pdf(pdwithout,x);
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Poisson''s Ratio')
title('Double-Mirror Symmetric Designs');
xlim([-0.6 0.6]);
set(gca,'fontsize', 12)

han=axes(fig6,'visible','off'); 
han.YLabel.Visible='on';
ylabel(han,'Probability Density');
set(gca,'fontsize', 12)

%%% SEVENTH PLOT: Normalized C66 vs Stacked Members
fig7 = figure(7);

asCval_with_6 = (2.*asC66_with_6)./(asC11_with_6+asC22_with_6);
asCval_without_6 = (2.*asC66_without_6)./(asC11_without_6+asC22_without_6);

fsCval_with_6 = (2.*fsC66_with_6)./(fsC11_with_6+fsC22_with_6);
fsCval_without_6 = (2.*fsC66_without_6)./(fsC11_without_6+fsC22_without_6);

sCval_without_6 = (2.*sC66Basket)./(sC11Basket+sC22Basket);

% First plot
subplot(3,1,1); 
x = 0:0.001:0.3;
[~,~,h_with] = ksdensity((asCval_with_6));
pdwith = fitdist((asCval_with_6),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity((asCval_without_6));
pdwithout = fitdist((asCval_without_6),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
title('Asymmetric Designs');
xlabel('Normalized Shear Stiffness')
legend('With Stacks','Without Stacks');
xlim([0 0.3]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Second plot
subplot(3,1,2); 
x = 0:0.001:0.3;
[~,~,h_with] = ksdensity((fsCval_with_6));
pdwith = fitdist((fsCval_with_6),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity((fsCval_without_6));
pdwithout = fitdist((fsCval_without_6),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normalized Shear Stiffness')
title('Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Third plot
subplot(3,1,3); 
x = 0:0.001:0.3;
[~,~,h_without] = ksdensity((sCval_without_6));
pdwithout = fitdist((sCval_without_6)','Kernel','Width',h_without);
ywithout = pdf(pdwithout,x);
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normalized Shear Stiffness')
title('Double-Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'fontsize', 12)

han=axes(fig7,'visible','off'); 
han.YLabel.Visible='on';
ylabel(han,'Probability Density');
set(gca,'fontsize', 12)

%%% EIGTH PLOT: v21 vs Stacked Members
fig8 = figure(8);

% First plot
subplot(3,1,1); 
x = -1:0.001:1;
[~,~,h_with] = ksdensity(asvyx_with_6);
pdwith = fitdist(asvyx_with_6,'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(asvyx_without_6);
pdwithout = fitdist(asvyx_without_6,'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
title('Asymmetric Designs');
xlabel('Poisson''s Ratio')
legend('With Stacks','Without Stacks','Location','northwest');
xlim([-0.6 0.6]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Second plot
subplot(3,1,2); 
x = -1:0.001:1;
[~,~,h_with] = ksdensity(fsvyx_with_6);
pdwith = fitdist(fsvyx_with_6,'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(fsvyx_without_6);
pdwithout = fitdist(fsvyx_without_6,'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Poisson''s Ratio')
title('Mirror Symmetric Designs');
xlim([-0.6 0.6]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Third plot
subplot(3,1,3); 
x = -1:0.001:1;
[~,~,h_without] = ksdensity(svyxBasket);
pdwithout = fitdist(svyxBasket','Kernel','Width',h_without);
ywithout = pdf(pdwithout,x);
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Poisson''s Ratio')
title('Double-Mirror Symmetric Designs');
xlim([-0.6 0.6]);
set(gca,'fontsize', 12)

%sgtitle('Effect of Stacked Members on Poisson''s Ratio')
han=axes(fig8,'visible','off'); 
han.YLabel.Visible='on';
ylabel(han,'Probability Density');
set(gca,'fontsize', 12)

%% Kernel-Fit Plots of Design Characteristics vs Normal Stiffness (3x3)

%%% FIRST PLOT: Normal Stiffness vs Pivot Points
fig1 = figure(1);

% First plot
subplot(3,1,1); 
x = 0:0.001:1;
[~,~,h_with] = ksdensity(([asC11_with_2;asC22_with_2]./E));
pdwith = fitdist(([asC11_with_2;asC22_with_2]./E),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(([asC11_without_2;asC22_without_2]./E));
pdwithout = fitdist(([asC11_without_2;asC22_without_2]./E),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normal Stiffness')
legend('With Pivot Points','Without Pivot Points');
title('Asymmetric Designs');
xlim([0 0.3]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Second plot
subplot(3,1,2); 
x = 0:0.001:1;
[~,~,h_with] = ksdensity(([fsC11_with_2;fsC22_with_2]./E));
pdwith = fitdist(([fsC11_with_2;fsC22_with_2]./E),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(([fsC11_without_2;fsC22_without_2]./E));
pdwithout = fitdist(([fsC11_without_2;fsC22_without_2]./E),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normal Stiffness')
title('Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'fontsize', 12)

% Third plot
subplot(3,1,3); 
x = 0:0.001:1;
[~,~,h_with] = ksdensity(([sC11_with_2;sC22_with_2]./E));
pdwith = fitdist(([sC11_with_2;sC22_with_2]./E),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(([sC11_without_2;sC22_without_2]./E));
pdwithout = fitdist(([sC11_without_2;sC22_without_2]./E),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normal Stiffness')
title('Double-Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

han=axes(fig1,'visible','off'); 
han.YLabel.Visible='on';
ylabel(han,'Probability Density');
set(gca,'fontsize', 12)

%%% SECOND PLOT: Normal Stiffness vs Arrows
fig2 = figure(2);

% First plot
subplot(3,1,1); 
x = 0:0.001:1;
[~,~,h_with] = ksdensity(([asC11_with_4;asC22_with_4]./E));
pdwith = fitdist(([asC11_with_4;asC22_with_4]./E),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(([asC11_without_4;asC22_without_4]./E));
pdwithout = fitdist(([asC11_without_4;asC22_without_4]./E),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normal Stiffness')
legend('With Arrows','Without Arrows');
title('Asymmetric Designs');
xlim([0 0.3]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Second plot
subplot(3,1,2); 
x = 0:0.001:1;
[~,~,h_with] = ksdensity(([fsC11_with_4;fsC22_with_4]./E));
pdwith = fitdist(([fsC11_with_4;fsC22_with_4]./E),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(([fsC11_without_4;fsC22_without_4]./E));
pdwithout = fitdist(([fsC11_without_4;fsC22_without_4]./E),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normal Stiffness')
title('Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'fontsize', 12)

% Third plot
subplot(3,1,3); 
x = 0:0.001:1;
[~,~,h_without] = ksdensity(([sC11_without_4;sC22_without_4]./E));
pdwithout = fitdist(([sC11_without_4;sC22_without_4]./E),'Kernel','Width',h_without);
ywithout = pdf(pdwithout,x);
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normal Stiffness')
title('Double-Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'fontsize', 12)

han=axes(fig2,'visible','off'); 
han.YLabel.Visible='on';
ylabel(han,'Probability Density');
set(gca,'fontsize', 12)

%%% THIRD PLOT: Normal Stiffness vs Spider Nodes
fig3 = figure(3);

% First plot
subplot(3,1,1); 
x = 0:0.001:1;
[~,~,h_with] = ksdensity(([asC11_with_5;asC22_with_5]./E));
pdwith = fitdist(([asC11_with_5;asC22_with_5]./E),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(([asC11_without_5;asC22_without_5]./E));
pdwithout = fitdist(([asC11_without_5;asC22_without_5]./E),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normal Stiffness')
legend('With Spider Nodes','Without Spider Nodes');
title('Asymmetric Designs');
xlim([0 0.3]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Second plot
subplot(3,1,2); 
x = 0:0.001:1;
[~,~,h_with] = ksdensity(([fsC11_with_5;fsC22_with_5]./E));
pdwith = fitdist(([fsC11_with_5;fsC22_with_5]./E),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(([fsC11_without_5;fsC22_without_5]./E));
pdwithout = fitdist(([fsC11_without_5;fsC22_without_5]./E),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normal Stiffness')
title('Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'fontsize', 12)

% Third plot
subplot(3,1,3); 
x = 0:0.001:1;
[~,~,h_without] = ksdensity(([sC11_without_5;sC22_without_5]./E));
pdwithout = fitdist(([sC11_without_5;sC22_without_5]./E),'Kernel','Width',h_without);
ywithout = pdf(pdwithout,x);
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normal Stiffness')
title('Double-Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'fontsize', 12)

han=axes(fig3,'visible','off'); 
han.YLabel.Visible='on';
ylabel(han,'Probability Density');
set(gca,'fontsize', 12)

%%% FOURTH PLOT: Normal Stiffness vs Stacked Members
fig4 = figure(4);

% First plot
subplot(3,1,1); 
x = 0:0.001:1;
[~,~,h_with] = ksdensity(([asC11_with_6;asC22_with_6]./E));
pdwith = fitdist(([asC11_with_6;asC22_with_6]./E),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(([asC11_without_6;asC22_without_6]./E));
pdwithout = fitdist(([asC11_without_6;asC22_without_6]./E),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normal Stiffness')
legend('With Stacks','Without Stacks');
title('Asymmetric Designs');
xlim([0 0.3]);
set(gca,'fontsize', 12)

% Second plot
subplot(3,1,2); 
x = 0:0.001:1;
[~,~,h_with] = ksdensity(([fsC11_with_6;fsC22_with_6]./E));
pdwith = fitdist(([fsC11_with_6;fsC22_with_6]./E),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(([fsC11_without_6;fsC22_without_6]./E));
pdwithout = fitdist(([fsC11_without_6;fsC22_without_6]./E),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normal Stiffness')
title('Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'fontsize', 12)

% Third plot
subplot(3,1,3); 
x = 0:0.001:1;
[~,~,h_without] = ksdensity(([sC11_without_6;sC22_without_6]./E));
pdwithout = fitdist(([sC11_without_6;sC22_without_6]./E),'Kernel','Width',h_without);
ywithout = pdf(pdwithout,x);
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normal Stiffness')
title('Double-Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'fontsize', 12)

han=axes(fig4,'visible','off'); 
han.YLabel.Visible='on';
ylabel(han,'Probability Density');
set(gca,'fontsize', 12)

%% Two-Sample K-S Tests (3x3)

% Initialize values
KStestResults = zeros(3,12);

% First test: Pivot Points vs. C66
KStestResults(1,1) = kstest2(asCval_with_2,asCval_without_2);
KStestResults(2,1) = kstest2(fsCval_with_2,fsCval_without_2);
KStestResults(3,1) = kstest2(sCval_with_2,sCval_without_2);

% Second test: Pivot Points vs. v21
KStestResults(1,2) = kstest2(asvyx_with_2,asvyx_without_2);
KStestResults(2,2) = kstest2(fsvyx_with_2,fsvyx_without_2);
KStestResults(3,2) = kstest2(svyx_with_2,svyx_without_2);

% Third test: Pivot Points vs. Normal stiffness
KStestResults(1,3) = kstest2([asC11_with_2;asC22_with_2]./E,[asC11_without_2;asC22_without_2]./E);
KStestResults(2,3) = kstest2([fsC11_with_2;fsC22_with_2]./E,[fsC11_without_2;fsC22_without_2]./E);
KStestResults(3,3) = kstest2([sC11_with_2;sC22_with_2]./E,[sC11_without_2;sC22_without_2]./E);

% Fourth test: Arrows vs. C66
KStestResults(1,4) = kstest2(asCval_with_4,asCval_without_4);
KStestResults(2,4) = kstest2(fsCval_with_4,fsCval_without_4);

% Fifth test: Arrows vs. v21
KStestResults(1,5) = kstest2(asvyx_with_4,asvyx_without_4);
KStestResults(2,5) = kstest2(fsvyx_with_4,fsvyx_without_4);

% Sixth test: Arrows vs. Normal stiffness
KStestResults(1,6) = kstest2([asC11_with_4;asC22_with_4]./E,[asC11_without_6;asC22_without_4]./E);
KStestResults(2,6) = kstest2([fsC11_with_4;fsC22_with_4]./E,[fsC11_without_6;fsC22_without_4]./E);

% Seventh test: Spider Nodes vs. C66
KStestResults(1,7) = kstest2(asCval_with_5,asCval_without_5);
KStestResults(2,7) = kstest2(fsCval_with_5,fsCval_without_5);

% Eigth test: Spider Nodes vs. v21
KStestResults(1,8) = kstest2(asvyx_with_5,asvyx_without_5);
KStestResults(2,8) = kstest2(fsvyx_with_5,fsvyx_without_5);

% Ninth test: Spider Nodes vs. Normal stiffness
KStestResults(1,9) = kstest2([asC11_with_5;asC22_with_5]./E,[asC11_without_6;asC22_without_5]./E);
KStestResults(2,9) = kstest2([fsC11_with_5;fsC22_with_5]./E,[fsC11_without_6;fsC22_without_5]./E);

% Tenth test: Stacked Members vs. C66
KStestResults(1,10) = kstest2(asCval_with_6,asCval_without_6);
KStestResults(2,10) = kstest2(fsCval_with_6,fsCval_without_6);

% Eleventh test: Stacked Members vs. v21
KStestResults(1,11) = kstest2(asvyx_with_6,asvyx_without_6);
KStestResults(2,11) = kstest2(fsvyx_with_6,fsvyx_without_6);

% Twelfth test: Stacked Members vs. Normal stiffness
KStestResults(1,12) = kstest2([asC11_with_6;asC22_with_6]./E,[asC11_without_6;asC22_without_6]./E);
KStestResults(2,12) = kstest2([fsC11_with_6;fsC22_with_6]./E,[fsC11_without_6;fsC22_without_6]./E);

%% Kernel-Fit Plots of Design Characteristics vs Properties (5x5)

%%% FIRST PLOT: Normalized C66 vs Pivot Points
fig1 = figure(1);

asCval_with_2 = (2.*asC66_with_2)./(asC11_with_2+asC22_with_2);
asCval_without_2 = (2.*asC66_without_2)./(asC11_without_2+asC22_without_2);

fsCval_with_2 = (2.*fsC66_with_2)./(fsC11_with_2+fsC22_with_2);
fsCval_without_2 = (2.*fsC66_without_2)./(fsC11_without_2+fsC22_without_2);

sCval_with_2 = (2.*sC66_with_2)./(sC11_with_2+sC22_with_2);
sCval_without_2 = (2.*sC66_without_2)./(sC11_without_2+sC22_without_2);

% First plot
subplot(3,1,1); 
x = 0:0.001:0.3;
[~,~,h_with] = ksdensity((asCval_with_2));
pdwith = fitdist((asCval_with_2),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity((asCval_without_2));
pdwithout = fitdist((asCval_without_2),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normalized Shear Stiffness')
legend('With Pivot Points','Without Pivot Points');
title('Asymmetric Designs');
xlim([0 0.3]);
set(gca,'fontsize', 12)

% Second plot
subplot(3,1,2); 
x = 0:0.001:0.3;
[~,~,h_with] = ksdensity((fsCval_with_2));
pdwith = fitdist((fsCval_with_2),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity((fsCval_without_2));
pdwithout = fitdist((fsCval_without_2),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normalized Shear Stiffness')
title('Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'fontsize', 12)

% Third plot
subplot(3,1,3); 
x = 0:0.001:0.3;
[~,~,h_with] = ksdensity((sCval_with_2));
pdwith = fitdist((sCval_with_2),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity((sCval_without_2));
pdwithout = fitdist((sCval_without_2),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normalized Shear Stiffness')
title('Double-Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

han=axes(fig1,'visible','off'); 
han.YLabel.Visible='on';
ylabel(han,'Probability Density');
set(gca,'fontsize', 12)

%%% SECOND PLOT: v21 vs Pivot Points
fig2 = figure(2);

% First plot
subplot(3,1,1); 
x = -1:0.001:1;
[~,~,h_with] = ksdensity(asvyx_with_2);
pdwith = fitdist(asvyx_with_2,'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(asvyx_without_2);
pdwithout = fitdist(asvyx_without_2,'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Poisson''s Ratio')
legend('With Pivot Points','Without Pivot Points','Location','northwest');
title('Asymmetric Designs');
xlim([-0.6 0.6]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Second plot
subplot(3,1,2); 
x = -1:0.001:1;
[~,~,h_with] = ksdensity(fsvyx_with_2);
pdwith = fitdist(fsvyx_with_2,'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(fsvyx_without_2);
pdwithout = fitdist(fsvyx_without_2,'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
title('Mirror Symmetric Designs');
xlabel('Poisson''s Ratio')
xlim([-0.6 0.6]);
set(gca,'fontsize', 12)

% Third plot
subplot(3,1,3); 
x = -1:0.001:1;
[~,~,h_with] = ksdensity(svyx_with_2);
pdwith = fitdist(svyx_with_2,'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(svyx_without_2);
pdwithout = fitdist(svyx_without_2,'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
title('Double-Mirror Symmetric Designs');
xlabel('Poisson''s Ratio')
xlim([-0.6 0.6]);
set(gca,'fontsize', 12)

han=axes(fig2,'visible','off'); 
han.YLabel.Visible='on';
ylabel(han,'Probability Density');
set(gca,'fontsize', 12)

%%% THIRD PLOT: Normalized C66 vs Arrows
fig3 = figure(3);

asCval_with_4 = (2.*asC66_with_4)./(asC11_with_4+asC22_with_4);
asCval_without_4 = (2.*asC66_without_4)./(asC11_without_4+asC22_without_4);

fsCval_with_4 = (2.*fsC66_with_4)./(fsC11_with_4+fsC22_with_4);
fsCval_without_4 = (2.*fsC66_without_4)./(fsC11_without_4+fsC22_without_4);

sCval_with_4 = (2.*sC66_with_4)./(sC11_with_4+sC22_with_4);
sCval_without_4 = (2.*sC66_without_4)./(sC11_without_4+sC22_without_4);

% First plot
subplot(3,1,1); 
x = 0:0.001:0.3;
[~,~,h_with] = ksdensity((asCval_with_4));
pdwith = fitdist((asCval_with_4),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity((asCval_without_4));
pdwithout = fitdist((asCval_without_4),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
title('Asymmetric Designs');
xlabel('Normalized Shear Stiffness')
legend('With Arrows','Without Arrows');
xlim([0 0.3]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Second plot
subplot(3,1,2); 
x = 0:0.001:0.3;
[~,~,h_with] = ksdensity((fsCval_with_4));
pdwith = fitdist((fsCval_with_4),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity((fsCval_without_4));
pdwithout = fitdist((fsCval_without_4),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normalized Shear Stiffness')
title('Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Third plot
subplot(3,1,3); 
x = 0:0.001:0.3;
[~,~,h_with] = ksdensity((sCval_with_4));
pdwith = fitdist((sCval_with_4),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity((sCval_without_4));
pdwithout = fitdist((sCval_without_4),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normalized Shear Stiffness')
title('Double-Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

han=axes(fig3,'visible','off'); 
han.YLabel.Visible='on';
ylabel(han,'Probability Density');
set(gca,'fontsize', 12)

%%% FOURTH PLOT: v21 vs Arrows
fig4 = figure(4);

% First plot
subplot(3,1,1); 
x = -1:0.001:1;
[~,~,h_with] = ksdensity(asvyx_with_4);
pdwith = fitdist(asvyx_with_4,'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(asvyx_without_4);
pdwithout = fitdist(asvyx_without_4,'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
title('Asymmetric Designs');
xlabel('Poisson''s Ratio')
legend('With Arrows','Without Arrows','Location','northwest');
xlim([-0.6 0.6]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Second plot
subplot(3,1,2); 
x = -1:0.001:1;
[~,~,h_with] = ksdensity(fsvyx_with_4);
pdwith = fitdist(fsvyx_with_4,'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(fsvyx_without_4);
pdwithout = fitdist(fsvyx_without_4,'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Poisson''s Ratio')
title('Mirror Symmetric Designs');
xlim([-0.6 0.6]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Third plot
subplot(3,1,3); 
x = -1:0.001:1;
[~,~,h_with] = ksdensity(svyx_with_4);
pdwith = fitdist(svyx_with_4,'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(svyx_without_4);
pdwithout = fitdist(svyx_without_4,'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Poisson''s Ratio')
title('Double-Mirror Symmetric Designs');
xlim([-0.6 0.6]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

han=axes(fig4,'visible','off'); 
han.YLabel.Visible='on';
ylabel(han,'Probability Density');
set(gca,'fontsize', 12)

%%% FIFTH PLOT: Normalized C66 vs Spider Nodes
fig5 = figure(5);

asCval_with_5 = (2.*asC66_with_5)./(asC11_with_5+asC22_with_5);
asCval_without_5 = (2.*asC66_without_5)./(asC11_without_5+asC22_without_5);

fsCval_with_5 = (2.*fsC66_with_5)./(fsC11_with_5+fsC22_with_5);
fsCval_without_5 = (2.*fsC66_without_5)./(fsC11_without_5+fsC22_without_5);

sCval_with_5 = (2.*sC66_with_5)./(sC11_with_5+sC22_with_5);
sCval_without_5 = (2.*sC66_without_5)./(sC11_without_5+sC22_without_5);

% First plot
subplot(3,1,1); 
x = 0:0.001:0.3;
[~,~,h_with] = ksdensity((asCval_with_5));
pdwith = fitdist((asCval_with_5),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity((asCval_without_5));
pdwithout = fitdist((asCval_without_5),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
title('Asymmetric Designs');
xlabel('Normalized Shear Stiffness')
legend('With Spider Nodes','Without Spider Nodes');
xlim([0 0.3]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Second plot
subplot(3,1,2); 
x = 0:0.001:0.3;
[~,~,h_with] = ksdensity((fsCval_with_5));
pdwith = fitdist((fsCval_with_5),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity((fsCval_without_5));
pdwithout = fitdist((fsCval_without_5),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normalized Shear Stiffness')
title('Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Third plot
subplot(3,1,3); 
x = 0:0.001:0.3;
[~,~,h_with] = ksdensity((sCval_with_5));
pdwith = fitdist((sCval_with_5),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity((sCval_without_5));
pdwithout = fitdist((sCval_without_5),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normalized Shear Stiffness')
title('Double-Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

han=axes(fig5,'visible','off'); 
han.YLabel.Visible='on';
ylabel(han,'Probability Density');
set(gca,'fontsize', 12)



%%% SIXTH PLOT: v21 vs Spider Nodes
fig6 = figure(6);

% First plot
subplot(3,1,1); 
x = -1:0.001:1;
[~,~,h_with] = ksdensity(asvyx_with_5);
pdwith = fitdist(asvyx_with_5,'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(asvyx_without_5);
pdwithout = fitdist(asvyx_without_5,'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
title('Asymmetric Designs');
xlabel('Poisson''s Ratio')
legend('With Spider Nodes','Without Spider Nodes','Location','northwest');
xlim([-0.6 0.6]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Second plot
subplot(3,1,2); 
x = -1:0.001:1;
[~,~,h_with] = ksdensity(fsvyx_with_5);
pdwith = fitdist(fsvyx_with_5,'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(fsvyx_without_5);
pdwithout = fitdist(fsvyx_without_5,'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Poisson''s Ratio')
title('Mirror Symmetric Designs');
xlim([-0.6 0.6]);
set(gca,'fontsize', 12)

% Third plot
subplot(3,1,3); 
x = -1:0.001:1;
[~,~,h_with] = ksdensity(svyx_with_5);
pdwith = fitdist(svyx_with_5,'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(svyx_without_5);
pdwithout = fitdist(svyx_without_5,'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Poisson''s Ratio')
title('Double-Mirror Symmetric Designs');
xlim([-0.6 0.6]);
set(gca,'fontsize', 12)

han=axes(fig6,'visible','off'); 
han.YLabel.Visible='on';
ylabel(han,'Probability Density');
set(gca,'fontsize', 12)

%%% SEVENTH PLOT: Normalized C66 vs Stacked Members
fig7 = figure(7);

asCval_with_6 = (2.*asC66_with_6)./(asC11_with_6+asC22_with_6);
asCval_without_6 = (2.*asC66_without_6)./(asC11_without_6+asC22_without_6);

fsCval_with_6 = (2.*fsC66_with_6)./(fsC11_with_6+fsC22_with_6);
fsCval_without_6 = (2.*fsC66_without_6)./(fsC11_without_6+fsC22_without_6);

sCval_with_6 = (2.*sC66_with_6)./(sC11_with_6+sC22_with_6);
sCval_without_6 = (2.*sC66_without_6)./(sC11_without_6+sC22_without_6);

% First plot
subplot(3,1,1); 
x = 0:0.001:0.3;
[~,~,h_with] = ksdensity((asCval_with_6));
pdwith = fitdist((asCval_with_6),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity((asCval_without_6));
pdwithout = fitdist((asCval_without_6),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
title('Asymmetric Designs');
xlabel('Normalized Shear Stiffness')
legend('With Stacks','Without Stacks');
xlim([0 0.3]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Second plot
subplot(3,1,2); 
x = 0:0.001:0.3;
[~,~,h_with] = ksdensity((fsCval_with_6));
pdwith = fitdist((fsCval_with_6),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity((fsCval_without_6));
pdwithout = fitdist((fsCval_without_6),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normalized Shear Stiffness')
title('Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Third plot
subplot(3,1,3); 
x = 0:0.001:0.3;
[~,~,h_with] = ksdensity((sCval_with_6));
pdwith = fitdist((sCval_with_6),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity((sCval_without_6));
pdwithout = fitdist((sCval_without_6),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normalized Shear Stiffness')
title('Double-Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

han=axes(fig7,'visible','off'); 
han.YLabel.Visible='on';
ylabel(han,'Probability Density');
set(gca,'fontsize', 12)

%%% EIGTH PLOT: v21 vs Stacked Members
fig8 = figure(8);

% First plot
subplot(3,1,1); 
x = -1:0.001:1;
[~,~,h_with] = ksdensity(asvyx_with_6);
pdwith = fitdist(asvyx_with_6,'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(asvyx_without_6);
pdwithout = fitdist(asvyx_without_6,'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
title('Asymmetric Designs');
xlabel('Poisson''s Ratio')
legend('With Stacks','Without Stacks','Location','northwest');
xlim([-0.6 0.6]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Second plot
subplot(3,1,2); 
x = -1:0.001:1;
[~,~,h_with] = ksdensity(fsvyx_with_6);
pdwith = fitdist(fsvyx_with_6,'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(fsvyx_without_6);
pdwithout = fitdist(fsvyx_without_6,'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Poisson''s Ratio')
title('Mirror Symmetric Designs');
xlim([-0.6 0.6]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Third plot
subplot(3,1,3); 
x = -1:0.001:1;
[~,~,h_with] = ksdensity(svyx_with_6);
pdwith = fitdist(svyx_with_6,'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(svyx_without_6);
pdwithout = fitdist(svyx_without_6,'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Poisson''s Ratio')
title('Double-Mirror Symmetric Designs');
xlim([-0.6 0.6]);
set(gca,'fontsize', 12)

han=axes(fig8,'visible','off'); 
han.YLabel.Visible='on';
ylabel(han,'Probability Density');
set(gca,'fontsize', 12)

%% Kernel-Fit Plots of Design Characteristics vs Normal Stiffness (5x5)

%%% FIRST PLOT: Normal Stiffness vs Pivot Points
fig1 = figure(1);

% First plot
subplot(3,1,1); 
x = 0:0.001:1;
[~,~,h_with] = ksdensity(([asC11_with_2;asC22_with_2]./E));
pdwith = fitdist(([asC11_with_2;asC22_with_2]./E),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(([asC11_without_2;asC22_without_2]./E));
pdwithout = fitdist(([asC11_without_2;asC22_without_2]./E),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normal Stiffness')
legend('With Pivot Points','Without Pivot Points');
title('Asymmetric Designs');
xlim([0 0.3]);
set(gca,'fontsize', 12)

% Second plot
subplot(3,1,2); 
x = 0:0.001:1;
[~,~,h_with] = ksdensity(([fsC11_with_2;fsC22_with_2]./E));
pdwith = fitdist(([fsC11_with_2;fsC22_with_2]./E),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(([fsC11_without_2;fsC22_without_2]./E));
pdwithout = fitdist(([fsC11_without_2;fsC22_without_2]./E),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normal Stiffness')
title('Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'fontsize', 12)

% Third plot
subplot(3,1,3); 
x = 0:0.001:1;
[~,~,h_with] = ksdensity(([sC11_with_2;sC22_with_2]./E));
pdwith = fitdist(([sC11_with_2;sC22_with_2]./E),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(([sC11_without_2;sC22_without_2]./E));
pdwithout = fitdist(([sC11_without_2;sC22_without_2]./E),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normal Stiffness')
title('Double-Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'fontsize', 12)

han=axes(fig1,'visible','off'); 
han.YLabel.Visible='on';
ylabel(han,'Probability Density');
set(gca,'fontsize', 12)

%%% SECOND PLOT: Normal Stiffness vs Arrows
fig2 = figure(2);

% First plot
subplot(3,1,1); 
x = 0:0.001:1;
[~,~,h_with] = ksdensity(([asC11_with_4;asC22_with_4]./E));
pdwith = fitdist(([asC11_with_4;asC22_with_4]./E),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(([asC11_without_4;asC22_without_4]./E));
pdwithout = fitdist(([asC11_without_4;asC22_without_4]./E),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normal Stiffness')
legend('With Arrows','Without Arrows');
title('Asymmetric Designs');
xlim([0 0.3]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Second plot
subplot(3,1,2); 
x = 0:0.001:1;
[~,~,h_with] = ksdensity(([fsC11_with_4;fsC22_with_4]./E));
pdwith = fitdist(([fsC11_with_4;fsC22_with_4]./E),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(([fsC11_without_4;fsC22_without_4]./E));
pdwithout = fitdist(([fsC11_without_4;fsC22_without_4]./E),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normal Stiffness')
title('Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Third plot
subplot(3,1,3); 
x = 0:0.001:1;
[~,~,h_with] = ksdensity(([sC11_with_4;sC22_with_4]./E));
pdwith = fitdist(([sC11_with_4;sC22_with_4]./E),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(([sC11_without_4;sC22_without_4]./E));
pdwithout = fitdist(([sC11_without_4;sC22_without_4]./E),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normal Stiffness')
title('Double-Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

han=axes(fig2,'visible','off'); 
han.YLabel.Visible='on';
ylabel(han,'Probability Density');
set(gca,'fontsize', 12)

%%% THIRD PLOT: Normal Stiffness vs Spider Nodes
fig3 = figure(3);

% First plot
subplot(3,1,1); 
x = 0:0.001:1;
[~,~,h_with] = ksdensity(([asC11_with_5;asC22_with_5]./E));
pdwith = fitdist(([asC11_with_5;asC22_with_5]./E),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(([asC11_without_5;asC22_without_5]./E));
pdwithout = fitdist(([asC11_without_5;asC22_without_5]./E),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normal Stiffness')
legend('With Spider Nodes','Without Spider Nodes');
title('Asymmetric Designs');
xlim([0 0.3]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Second plot
subplot(3,1,2); 
x = 0:0.001:1;
[~,~,h_with] = ksdensity(([fsC11_with_5;fsC22_with_5]./E));
pdwith = fitdist(([fsC11_with_5;fsC22_with_5]./E),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(([fsC11_without_5;fsC22_without_5]./E));
pdwithout = fitdist(([fsC11_without_5;fsC22_without_5]./E),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normal Stiffness')
title('Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Third plot
subplot(3,1,3); 
x = 0:0.001:1;
[~,~,h_with] = ksdensity(([sC11_with_5;sC22_with_5]./E));
pdwith = fitdist(([sC11_with_5;sC22_with_5]./E),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(([sC11_without_5;sC22_without_5]./E));
pdwithout = fitdist(([sC11_without_5;sC22_without_5]./E),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normal Stiffness')
title('Double-Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

han=axes(fig3,'visible','off'); 
han.YLabel.Visible='on';
ylabel(han,'Probability Density');
set(gca,'fontsize', 12)

%%% FOURTH PLOT: Normal Stiffness vs Stacked Members
fig4 = figure(4);

% First plot
subplot(3,1,1); 
x = 0:0.001:1;
[~,~,h_with] = ksdensity(([asC11_with_6;asC22_with_6]./E));
pdwith = fitdist(([asC11_with_6;asC22_with_6]./E),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(([asC11_without_6;asC22_without_6]./E));
pdwithout = fitdist(([asC11_without_6;asC22_without_6]./E),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normal Stiffness')
legend('With Stacks','Without Stacks');
title('Asymmetric Designs');
xlim([0 0.3]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Second plot
subplot(3,1,2); 
x = 0:0.001:1;
[~,~,h_with] = ksdensity(([fsC11_with_6;fsC22_with_6]./E));
pdwith = fitdist(([fsC11_with_6;fsC22_with_6]./E),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(([fsC11_without_6;fsC22_without_6]./E));
pdwithout = fitdist(([fsC11_without_6;fsC22_without_6]./E),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normal Stiffness')
title('Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'Color',[0.8 0.975 0.9]);
set(gca,'fontsize', 12)

% Third plot
subplot(3,1,3); 
x = 0:0.001:1;
[~,~,h_with] = ksdensity(([sC11_with_6;sC22_with_6]./E));
pdwith = fitdist(([sC11_with_6;sC22_with_6]./E),'Kernel','Width',h_with);
[~,~,h_without] = ksdensity(([sC11_without_6;sC22_without_6]./E));
pdwithout = fitdist(([sC11_without_6;sC22_without_6]./E),'Kernel','Width',h_without);
ywith = pdf(pdwith,x);
ywithout = pdf(pdwithout,x);
plot(x,ywith,'k-','LineWidth',2)
hold on
plot(x,ywithout,'b-','LineWidth',2)
xlabel('Normal Stiffness')
title('Double-Mirror Symmetric Designs');
xlim([0 0.3]);
set(gca,'fontsize', 12)

han=axes(fig4,'visible','off'); 
han.YLabel.Visible='on';
ylabel(han,'Probability Density');
set(gca,'fontsize', 12)

%% Two-Sample K-S Tests (5x5)

% Initialize values
KStestResults = zeros(3,12);

% First test: Pivot Points vs. C66
KStestResults(1,1) = kstest2(asCval_with_2,asCval_without_2);
KStestResults(2,1) = kstest2(fsCval_with_2,fsCval_without_2);
KStestResults(3,1) = kstest2(sCval_with_2,sCval_without_2);

% Second test: Pivot Points vs. v21
KStestResults(1,2) = kstest2(asvyx_with_2,asvyx_without_2);
KStestResults(2,2) = kstest2(fsvyx_with_2,fsvyx_without_2);
KStestResults(3,2) = kstest2(svyx_with_2,svyx_without_2);

% Third test: Pivot Points vs. Normal Stiffness
KStestResults(1,3) = kstest2([asC11_with_2;asC22_with_2]./E,[asC11_without_2;asC22_without_2]./E);
KStestResults(2,3) = kstest2([fsC11_with_2;fsC22_with_2]./E,[fsC11_without_2;fsC22_without_2]./E);
KStestResults(3,3) = kstest2([sC11_with_2;sC22_with_2]./E,[sC11_without_2;sC22_without_2]./E);

% Fourth test: Arrows vs. C66
KStestResults(1,4) = kstest2(asCval_with_4,asCval_without_4);
KStestResults(2,4) = kstest2(fsCval_with_4,fsCval_without_4);
KStestResults(3,4) = kstest2(sCval_with_4,sCval_without_4);

% Fifth test: Arrows vs. v21
KStestResults(1,5) = kstest2(asvyx_with_4,asvyx_without_4);
KStestResults(2,5) = kstest2(fsvyx_with_4,fsvyx_without_4);
KStestResults(3,5) = kstest2(svyx_with_4,svyx_without_4);

% Sixth test: Arrows vs. Normal Stiffness
KStestResults(1,6) = kstest2([asC11_with_4;asC22_with_4]./E,[asC11_without_4;asC22_without_4]./E);
KStestResults(2,6) = kstest2([fsC11_with_4;fsC22_with_4]./E,[fsC11_without_4;fsC22_without_4]./E);
KStestResults(3,6) = kstest2([sC11_with_4;sC22_with_4]./E,[sC11_without_4;sC22_without_4]./E);

% Seventh test: Spider Nodes vs. C66
KStestResults(1,7) = kstest2(asCval_with_5,asCval_without_5);
KStestResults(2,7) = kstest2(fsCval_with_5,fsCval_without_5);
KStestResults(3,7) = kstest2(sCval_with_5,sCval_without_5);

% Eigth test: Spider Nodes vs. v21
KStestResults(1,8) = kstest2(asvyx_with_5,asvyx_without_5);
KStestResults(2,8) = kstest2(fsvyx_with_5,fsvyx_without_5);
KStestResults(3,8) = kstest2(svyx_with_5,svyx_without_5);

% Ninth test: Arrows vs. Normal Stiffness
KStestResults(1,9) = kstest2([asC11_with_5;asC22_with_5]./E,[asC11_without_5;asC22_without_5]./E);
KStestResults(2,9) = kstest2([fsC11_with_5;fsC22_with_5]./E,[fsC11_without_5;fsC22_without_5]./E);
KStestResults(3,9) = kstest2([sC11_with_5;sC22_with_5]./E,[sC11_without_5;sC22_without_5]./E);

% Tenth test: Stacked Members vs. C66
KStestResults(1,10) = kstest2(asCval_with_6,asCval_without_6);
KStestResults(2,10) = kstest2(fsCval_with_6,fsCval_without_6);
KStestResults(3,10) = kstest2(sCval_with_6,sCval_without_6);

% Eleventh test: Stacked Members vs. v21
KStestResults(1,11) = kstest2(asvyx_with_6,asvyx_without_6);
KStestResults(2,11) = kstest2(fsvyx_with_6,fsvyx_without_6);
KStestResults(3,11) = kstest2(svyx_with_6,svyx_without_6);

% Sixth test: Arrows vs. Normal Stiffness
KStestResults(1,12) = kstest2([asC11_with_6;asC22_with_6]./E,[asC11_without_6;asC22_without_6]./E);
KStestResults(2,12) = kstest2([fsC11_with_6;fsC22_with_6]./E,[fsC11_without_6;fsC22_without_6]./E);
KStestResults(3,12) = kstest2([sC11_with_6;sC22_with_6]./E,[sC11_without_6;sC22_without_6]./E);

%% ---------- %%
%%% SUBFUNCTIONS
% FUNCTION TO GENERATE NODAL COORDINATES BASED ON GRID SIZE
function NC = generateNC(sel,sidenum)
    notchvec = linspace(0,1,sidenum);
    NC = [];
    for i = 1:1:sidenum
        for j = 1:1:sidenum
            NC = [NC;notchvec(i),notchvec(j)];
        end
    end
    NC = sel.*NC;
end

% FUNCTION TO INPUT SYMMETRIC DESIGNS
function symmCA = inputSymmCA()
    symmCA(1) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4]};
    symmCA(2) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;2,5;5,8;4,5;5,6]};
    symmCA(3) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;1,5;3,5;5,7;5,9]};
    symmCA(4) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;2,4;2,6;6,8;4,8]};
    symmCA(5) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;2,5;5,8;4,5;5,6;1,5;3,5;5,7;5,9]};
    symmCA(6) = {[1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;2,5;5,8;4,5;5,6;2,4;2,6;6,8;4,8]};
    symmCA(7) = {[2,5;5,8;4,5;5,6]};
    symmCA(8) = {[1,5;3,5;5,7;5,9]};
    symmCA(9) = {[2,4;2,6;6,8;4,8]};
    symmCA(10) = {[2,5;5,8;4,5;5,6;1,5;3,5;5,7;5,9]};
    symmCA(11) = {[2,5;5,8;4,5;5,6;2,4;2,6;6,8;4,8]};
    symmCA(12) = {[1,2;3,6;8,9;4,7;2,5;5,8;4,5;5,6]};
    symmCA(13) = {[1,2;3,6;8,9;4,7;1,5;3,5;5,7;5,9]};
    symmCA(14) = {[1,2;3,6;8,9;4,7;2,4;2,6;6,8;4,8]};
    symmCA(15) = {[1,2;3,6;8,9;4,7;2,5;5,8;4,5;5,6;1,5;3,5;5,7;5,9]};
    symmCA(16) = {[1,2;3,6;8,9;4,7;2,5;5,8;4,5;5,6;2,4;2,6;6,8;4,8]};
    symmCA(17) = {[2,3;6,9;7,8;1,4;2,5;5,8;4,5;5,6]};
    symmCA(18) = {[2,3;6,9;7,8;1,4;1,5;3,5;5,7;5,9]};
    symmCA(19) = {[2,3;6,9;7,8;1,4;2,4;2,6;6,8;4,8]};
    symmCA(20) = {[2,3;6,9;7,8;1,4;2,5;5,8;4,5;5,6;1,5;3,5;5,7;5,9]};
    symmCA(21) = {[2,3;6,9;7,8;1,4;2,5;5,8;4,5;5,6;2,4;2,6;6,8;4,8]};
end

% FUNCTION TO CALCULATE DESIGN FEASIBILITY 
function feasibilityScore = feasChecker(NC,CA_des)     
    feasibilityScore = 1;

    % FIRST CONSTRAINT: members only intersect at nodes (no crossing)
    % Sort points from left to right by x-position
    SortedCA = sortrows(CA_des);
    
    % Develop 4xM matrix of line segment endpoint coordinates, where M is 
    %   the number of truss members.  Each row of format (x1,y1,x2,y2),
    %   where point 1 is leftmost, point 2 is rightmost
    PosA = [NC(SortedCA(:,1),1),NC(SortedCA(:,1),2),...
            NC(SortedCA(:,2),1),NC(SortedCA(:,2),2)];
    
    % Loop through each pair of elements
    for i = 1:1:size(PosA,1)
        for j = 1:1:size(PosA,1)
            % Determine whether the given pair of elements intersects
            intersect = findLineSegIntersection([PosA(i,1),PosA(i,2)],...
                        [PosA(i,3),PosA(i,4)],[PosA(j,1),PosA(j,2)],...
                        [PosA(j,3),PosA(j,4)]);
            
            % Throw an error, given an intersection
            if intersect == true
                feasibilityScore = feasibilityScore - 0.1;
                if feasibilityScore < 0.1
                    return
                end
            end
        end
    end
    
    % SECOND CONSTRAINT: Elements (of either the same or different lengths)
    %   cannot overlap
    % Loop through each element
    for k = 1:1:size(SortedCA,1)
        % Loop through each element again, to consider each possible pair 
        %   of elements
        for q = 1:1:size(SortedCA,1)
            % Check if both elements share a common startpoint
            if (NC(SortedCA(k,1),1) == NC(SortedCA(q,1),1)) && ...
                (NC(SortedCA(k,1),2) == NC(SortedCA(q,1),2))
                % Check if both elements have the same slope (and reject 
                %    the design if so)
                mk = (NC(SortedCA(k,2),2)-NC(SortedCA(k,1),2))/...
                     (NC(SortedCA(k,2),1)-NC(SortedCA(k,1),1));
                mq = (NC(SortedCA(q,2),2)-NC(SortedCA(q,1),2))/...
                     (NC(SortedCA(q,2),1)-NC(SortedCA(q,1),1));
                % If the same element is being compared twice, move on
                if k == q
                    continue
                elseif mk == mq
                   feasibilityScore = feasibilityScore - 0.1;
                   if feasibilityScore < 0.1
                       return
                   end
                end
            % Check if both elements share a common endpoint    
            elseif (NC(SortedCA(k,2),1) == NC(SortedCA(q,2),1)) && ...
                   (NC(SortedCA(k,2),2) == NC(SortedCA(q,2),2))    
                % Check if both elements have the same slope (and reject 
                %    the design if so)
                mk = (NC(SortedCA(k,2),2)-NC(SortedCA(k,1),2))/...
                     (NC(SortedCA(k,2),1)-NC(SortedCA(k,1),1));
                mq = (NC(SortedCA(q,2),2)-NC(SortedCA(q,1),2))/...
                     (NC(SortedCA(q,2),1)-NC(SortedCA(q,1),1));
                % If the same element is being compared twice, move on
                if k == q
                    continue
                elseif mk == mq
                   feasibilityScore = feasibilityScore - 0.1;
                   if feasibilityScore < 0.1
                       return
                   end
                end
            end
        end
    end
end

% FUNCTION TO DETERMINE PRESENCE OF INTERSECTION (FOR FEAS CONSTRAINT #1)
function intersect = findLineSegIntersection(p1,q1,p2,q2)
    if (findOrientation(p1,q1,p2) ~= findOrientation(p1,q1,q2))&&...
            (findOrientation(p2,q2,p1) ~= findOrientation(p2,q2,q1)) 
        if isequal(p1,p2) || isequal(q1,q2) || ...
                isequal(p1,q2) || isequal(q1,p2)
            % Intersection due to shared start or end point
            intersect = false; 
        else
            intersect = true;
        end
    else
        % Intersection not present
        intersect = false;
    end
end

% FUNCTION TO CALCULATE ORIENTATION FROM 3 POINTS (FOR FEAS CONSTRAINT #1)
function orientation = findOrientation(p,q,r)
    val = ((q(2)-p(2))*(r(1)-q(1)))-((q(1)-p(1))*(r(2)-q(2)));
    if val == 0
        orientation = 0;
    elseif val > 0
        orientation = 1;
    else
        orientation = 2;
    end
end

% FUNCTION TO CHECK REPEATABILITY
function repeatabilityBool = repChecker_2D_V1(CA,sidenum)
    % Checking left & right sides for match
    leftside = [];
    for i = 1:1:(sidenum-1)
        elem = [i,i+1];
        leftside = [leftside,ismember(elem,CA,'row')];
    end
    rightside = [];
    for i = ((sidenum^2)-sidenum+1):1:((sidenum^2)-1)
        elem = [i,i+1];
        rightside = [rightside,ismember(elem,CA,'row')];
    end
    if leftside == rightside
        repeatabilityBool = 1;
    else
        repeatabilityBool = 0;
        return
    end
    
    % Checking top & bottom sides for match
    topside = [];
    for i = sidenum:sidenum:((sidenum^2)-sidenum)
        elem = [i,i+sidenum];
        topside = [topside,ismember(elem,CA,'row')];
    end
    bottomside = [];
    for i = 1:sidenum:((sidenum^2)-(2*sidenum)+1)
        elem = [i,i+sidenum];
        bottomside = [bottomside,ismember(elem,CA,'row')];
    end
    if topside == bottomside
        repeatabilityBool = 1;
    else
        repeatabilityBool = 0;
    end
end

% FUNCTION TO CALCULATE EDGE-RELATED INPUTS
function [CA,edgelog] = edgeVals(CA,sidenum)
       
    % Identify four sets of edge nodes
    edge1nodes = 1:1:sidenum;
    edge2nodes = ((sidenum^2)-sidenum+1):1:(sidenum^2);
    edge3nodes = 1:sidenum:((sidenum^2)-(sidenum)+1);
    edge4nodes = sidenum:sidenum:(sidenum^2);
             
    % Identify edge members and edge types
    edge1connA = ismember(CA(:,1),edge1nodes);
    edge1connB = ismember(CA(:,2),edge1nodes);
    edge1log = (edge1connA & edge1connB); 
    edge1type = (edge1connA & edge1connB); 
    edge2connA = ismember(CA(:,1),edge2nodes);
    edge2connB = ismember(CA(:,2),edge2nodes);
    edge2log = (edge2connA & edge2connB); 
    edge2type = 2.*(edge2connA & edge2connB);
    edge3connA = ismember(CA(:,1),edge3nodes);
    edge3connB = ismember(CA(:,2),edge3nodes);
    edge3log = (edge3connA & edge3connB); 
    edge3type = 3.*(edge3connA & edge3connB);
    edge4connA = ismember(CA(:,1),edge4nodes);
    edge4connB = ismember(CA(:,2),edge4nodes);
    edge4log = (edge4connA & edge4connB); 
    edge4type = 4.*(edge4connA & edge4connB);
    
    edgelog = edge1log + edge2log + edge3log + edge4log;
    edgetype = edge1type + edge2type + edge3type + edge4type;
    
    % Modify order of edge members depending on type
    for i = 1:1:size(CA,1)
        if (edgetype(i) == 1) || (edgetype(i) == 4)
            if CA(i,2) > CA(i,1)
                newrow = [CA(i,2),CA(i,1)];
                CA(i,:) = newrow;
            end
        elseif (edgetype(i) == 2) || (edgetype(i) == 3)
            if CA(i,1) > CA(i,2)
                newrow = [CA(i,2),CA(i,1)];
                CA(i,:) = newrow;
            end
        end
    end
end

% FUNCTION TO CALCULATE CROSS-SECTIONAL PROPERTIES OF EDGE MEMBERS
function [csPropsA,csPropsB] = calcCSProps(r)
    A = 0.5*pi*(r^2); % Area of section
    Iyy = (pi*(r^4))/8; % Moment of inertia about the y axis
    Izz = (pi*(r^4))/8; % Moment of inertia about the z axis
    Iw = 0; % Warping constant
    J = (pi*(r^4))/4; % Torsional constant
    CGy = 0; % y coordinate of centroid
    CGz = (4*r)/(3*pi); % z coordinate of centroid
    Iyz = A*CGy*CGz; % Product of inertia
    SHy = 0; % y coordinate of shear center
    SHz = (4*r)/(3*pi); % z coordinate of shear center
    TKz = r; % Thickness along Z axis (maximum height)
    TKy = r/2; % Thickness along Y axis (maximum width)
    csPropsA = [A,Iyz,Iw,CGy,CGz,SHy,SHz,TKz,TKy];
    csPropsB = [Iyy,Izz,J];
end

% FUNCTION TO DETERMINE NODAL PAIRS FOR CONSTRAINTS EQUATIONS
function [NP,numLRCE,numTBCE] = findNodalPairs(sidenum,CA,NC)
    % Add up counters based on nodal connectivities 
    edgein = 0.5:1:(0.5+size(NC,1));
    [N,~] = histcounts(CA,edgein);
    
    % Determine vertical-edge nodal pairs
    NP = []; numLRCE = 0;
    for i = 1:1:sidenum
        if (N(i) == 0) && (N(i+(sidenum*(sidenum-1))) == 0)
            % Both nodes are unoccupied, no constraint equation required
        elseif (N(i) == 0) && (N(i+(sidenum*(sidenum-1))) ~= 0)
            % Left-edge node is unoccupied
            index = 0; counter = 0;
            while index == 0
                lower = 1 + counter; 
                upper = sidenum + counter;
                Nlocal = N(lower:upper);
                index = find(Nlocal,1,'first');
                counter = counter + sidenum;
            end
            NP = [NP;(index+(counter-sidenum)),(i+(sidenum*(sidenum-1)))];
            numLRCE = numLRCE + 1;
        elseif (N(i) ~= 0) && (N(i+(sidenum*(sidenum-1))) == 0)
            % Right-edge node is unoccupied
            index = 0; counter = 0;
            while index == 0
                lower = ((sidenum^2)-sidenum+1) - counter; 
                upper = (sidenum^2) - counter;
                Nlocal = N(lower:upper);
                index = find(Nlocal,1,'first');
                counter = counter + sidenum;
            end
            NP = [NP;i,(index+((sidenum^2)-sidenum)-(counter-sidenum))];
            numLRCE = numLRCE + 1;
        else
            % Both nodes are occupied, substitute nodes not required
            NP = [NP;i,(i+(sidenum*(sidenum-1)))];
            numLRCE = numLRCE + 1;
        end
    end
    
    % Determine horizontal-edge nodal pairs
    numTBCE = 0;
    for i = 1:sidenum:((sidenum^2)-sidenum+1)
        if (N(i) == 0) && (N(i+sidenum-1) == 0)
            % Both nodes are unoccupied, no constraint equation required
        elseif (N(i) == 0) && (N(i+sidenum-1) ~= 0)
            % Bottom-edge node is unoccupied
            index = 0; counter = 0;
            while index == 0
                lower = 1 + counter; 
                upper = (sidenum^2)-sidenum+1 + counter;
                Nlocal = N(lower:upper);
                index = find(Nlocal,1,'first');
                counter = counter + 1;
            end
            NP = [NP;((((index-1)*sidenum)+1)+(counter-1)),(i+sidenum-1)];
            numTBCE = numTBCE + 1;
        elseif (N(i) ~= 0) && (N(i+sidenum-1) == 0)
            % Top-edge node is unoccupied
            index = 0; counter = 0;
            while index == 0
                lower = sidenum - counter; 
                upper = (sidenum^2) - counter;
                Nlocal = N(lower:upper);
                index = find(Nlocal,1,'first');
                counter = counter + 1;
            end
            NP = [NP;i,((index*sidenum)-(counter-1))];
            numTBCE = numTBCE + 1;
        else
            % Both nodes are occupied, substitute nodes not required
            NP = [NP;i,(i+sidenum-1)];
            numTBCE = numTBCE + 1;
        end
    end
end

% FUNCTION TO CALCULATE EFFECTIVE STIFFNESS TENSOR, VOLUME FRACTION
function [C,volFrac] = calcCTruss(NC,csPropsA,csPropsB,edgelog,NP,...
                                  numLRCE,numTBCE,sidenum,sel,r,E,CA)
    % Iterate through once for each strain component
    for y = 1:1:3
        % Define strain vector: [e11, e22, e12]'
        strainvec = [0;0;0];

        % Set that component equal to a dummy value (0.01 strain), set all 
        % other values to zero
        strainvec(y) = 0.001; 
        strainvec(3) = strainvec(3)*2;
        
        % Populate Input Vector
        qsel = sel./(10^-3); qr = r./(10^-6); 
        qcsPropsA = csPropsA./(10^-6); qcsPropsB = csPropsB./(10^-15);
        input_vec=[qsel;qr;E;size(CA,1);y;sidenum;numLRCE;numTBCE;...
          qcsPropsA';qcsPropsB';edgelog;NP(:,1);NP(:,2);CA(:,1);CA(:,2)];

        % Define paths and input/output files 
        ANSYS_path=...
 'C:\Program Files\ANSYS Inc\v211\ansys\bin\winx64\ANSYS211';
        APDL_name='Truss2D_NxN.txt';
        input_file_name='para_in.txt';
        output_file_name='para_out.txt';
        ANSYS_path=strcat('"',ANSYS_path,'"');
        disp('Current stage: ');disp(y);

        % Write input vector to para_in.txt
        %writematrix(input_vec,input_file_name);
        dlmwrite(input_file_name,input_vec,'delimiter',' ',...
           'precision','%8.8f','newline','pc');
        % Call ANSYS APDL, perform FEA
        status = system(sprintf('%s -b -p aa_r -i %s -o out.txt',...
            ANSYS_path,APDL_name));
        disp('Simulation status: ');disp(status);
        % Read the results from para_out.txt
        output_vec=load(output_file_name)';
        F_x = output_vec(1);
        F_y = output_vec(2);
        F_xy = output_vec(3);
        FBank(y,:) = [F_x,F_y,F_xy];

        % Calculate stress vector
        stressvec = (1/(sel*2*r)).*[F_x;F_y;F_xy];

        % Use strain and stress vectors to solve for the corresponding row 
        % of the C matrix
        Cdummy = stressvec/strainvec;
        C(:,y) = Cdummy(:,y);
    end
    
    % Calculate volume fraction
    volFrac = calcVF_NxN_feasOnly(NC,CA,r,sel,sidenum);
end

% FUNCTION TO CALCULATE SYMMETRY SCORE
function symmScore = symmHeuristic_2D(CA,NC,sel)
    % Initialize score
    symmScore = 2;
    
    % Check for full symmetry
    for i = 1:1:size(CA,1)
        % Find the rotation of both endpoints by 90 degrees (shifted back
        % to the first quadrant)
        ry1 = round((NC(CA(i,1),1)),6); 
        ry2 = round((NC(CA(i,2),1)),6);
        rx1 = round((-NC(CA(i,1),2) + sel),6); 
        rx2 = round((-NC(CA(i,2),2) + sel),6);
        
        % Find nodal indices of reflected points
        Ax = find(NC(:,1)==rx1); Ay = find(NC(:,2)==ry1);
        A = intersect(Ax,Ay);
        Bx = find(NC(:,1)==rx2); By = find(NC(:,2)==ry2);
        B = intersect(Bx,By);
        
        % Determine whether reflected members exist in CA (and demerit
        % design if not)
        if (ismember([A,B],CA,'rows') == false) && ...
           (ismember([B,A],CA,'rows') == false)
            symmScore = 0;
            break;
        end
    end
    
    % Check for flip symmetry
    if symmScore == 0
        symmScore = 1; fsChecker = 2;
        
        % Vertical flip symmetry
        for i = 1:1:size(CA,1)
            % Find the flip of both endpoints about the x axis (shifted 
            % back to the first quadrant)
            rx1 = round((NC(CA(i,1),1)),6); 
            rx2 = round((NC(CA(i,2),1)),6);
            ry1 = round((-NC(CA(i,1),2) + sel),6); 
            ry2 = round((-NC(CA(i,2),2) + sel),6);

            % Find nodal indices of x-reflected points
            Fx = find(NC(:,1)==rx1); Fy = find(NC(:,2)==ry1);
            F = intersect(Fx,Fy);
            Gx = find(NC(:,1)==rx2); Gy = find(NC(:,2)==ry2);
            G = intersect(Gx,Gy);

            % Determine whether reflected members exist in CA (and demerit
            % design if not)
            if (ismember([F,G],CA,'rows') == false) && ...
               (ismember([G,F],CA,'rows') == false)
                fsChecker = fsChecker - 1;
                break;
            end
        end
        
        % Horizontal flip symmetry
        for i = 1:1:size(CA,1)
            % Find the flip of both endpoints about the x axis (shifted 
            % back to the first quadrant)
            ry1 = round((NC(CA(i,1),1)),6); 
            ry2 = round((NC(CA(i,2),1)),6);
            rx1 = round((-NC(CA(i,1),2) + sel),6); 
            rx2 = round((-NC(CA(i,2),2) + sel),6);

            % Find nodal indices of y-reflected points
            Jx = find(NC(:,1)==rx1); Jy = find(NC(:,2)==ry1);
            J = intersect(Jx,Jy);
            Dx = find(NC(:,1)==rx2); Dy = find(NC(:,2)==ry2);
            D = intersect(Dx,Dy);

            % Determine whether reflected members exist in CA (and demerit
            % design if not)
            if (ismember([J,D],CA,'rows') == false) && ...
               (ismember([D,J],CA,'rows') == false)
                fsChecker = fsChecker - 1;
                break;
            end
        end
        
        % Check for either flip symmetry
        if fsChecker == 0
            symmScore = 0;
        end
    end
end

