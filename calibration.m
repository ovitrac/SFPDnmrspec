% first and principal script to load spectrum data :
%   - dbpur/dbmix   database of susbtance (or mixture) NMR spectra (db sorted by substances) (see nmrloaddb or nmrloadascii)
%   - dbxpur/dbxmix database of susbtance (or mixture) NMR spectra (db sorted by fields) (see nmrloaddb or nmrloadascii)
%   - dbfit         database of substance NMR theroretical spectra (see nmrloaddbfit)
%   - dbcalibration database of calibration of substance NMR spectra
%   - dbcorrsub     database of correlation multiple of susbtance in dictionnary (see nmrloadmatrix)
%   - dbclass       database of classification of susbtance NMR spectra from pairwise correlation coefficient
%   - dbmask        database on mask of substances (see nmrloadmask)
% and save in .mat including all databases and nfo
% add calibration by mol of proton 
% 3/10/2014 : add db of ranking substances NMR spectra (dbclass) and db of multiple correlation of susbtance in dictionnary (dbcorrsub)            

%% load raw data and fitted data
close all, clear all
local = fullfile(find_path_toolbox('rmnspec'),'data_pur');
% [dbpur,dbxpur,dbmask,mask,dbmultiplet,dbmix,dbxmix] = nmrloadbase; % load data of susbtance spectrum and mixture spectrum
[dbmask,mask,dbmultiplet] = nmrloadmask('path',local); %load substance mask info
[dbpur,dbxpur] = nmrloadascii('pur'); % load substance spectrum
[dbmix,dbxmix] = nmrloadascii('path',fullfile(find_path_toolbox('rmnspec'),'data_mixture'),'odsfile','extract.ods','sheetname','extract','sheetnamenfo','nfoextract','primarykey','reference'); %,'method','nearest'); %load mixture spectrum
dbfitname = 'dbfit.mat';
dbfit = nmrloaddbfit('path',local,'dbname',dbfitname); % load dbfit
mol = fieldnames(rmfield(dbfit,'help')); nmol = length(mol);

%% calibration area=f(H)
conc = cellfun(@(m) dbpur.(m).concentration/dbpur.(m).Mw,mol); %conversion g/l to mol/l
dbcalibration = struct([]);
% pcalib = cell(nfitmol,1);
% Scalib = cell(nfitmol,1);
for i=1:nmol
    nROI = size(dbfit.(mol{i}).all,1);
    area = [dbfit.(mol{i}).all(:,1).area]';
    proton = dbpur.(mol{i}).H(1:nROI); 
    molproton = conc(i).*proton;
    [pproton,Sproton] = polyfit(proton,area,1);
    [pmolproton,Smolproton] = polyfit(molproton,area*dbpur.(mol{i}).normvalue,1);
    dbcalibration(1).(mol{i}).p = pproton;
    dbcalibration(1).(mol{i}).S = Sproton;
    dbcalibration(1).(mol{i}).pmolproton = pmolproton;
    dbcalibration(1).(mol{i}).Smolproton = Smolproton;
    dbcalibration(1).(mol{i}).area = area;
    dbcalibration(1).(mol{i}).H = proton;
    dbcalibration(1).(mol{i}).molH  = molproton;
    dbcalibration(1).(mol{i}).InH = dbpur.(mol{i}).In/(pproton(1)*sum(proton)) ;
end
% concatenate
tmp = cellfun(@(m) dbcalibration.(m).InH,mol,'UniformOutput',false);
dbxcalibration.InH = cat(2,tmp{:});

%% calibration common of all susbtances in databases date : 23012014
molcalib = fieldnames(rmfield(dbfit,{'help' 'PE' 'PP' 'PS' 'Acetophenone' 'HBP4' 'DPGME' 'MBOCA' 'MBOCAR1' 'MBOCAR4' 'MIT' 'DBP' 'BBP' 'Tributylacetylcitrate'}));
nmolcalib = length(molcalib);
areacalib = cell(nmolcalib,1);
molH = cell(nmolcalib,1);
r2fit = cell(nmolcalib,1);
for i = 1:nmolcalib
    areacalib{i} = dbcalibration.(molcalib{i}).area*dbpur.(molcalib{i}).normvalue;
    molH{i} = dbcalibration.(molcalib{i}).molH;
    nROI = size(dbfit.(molcalib{i}).all,1); 
    r2fit{i} =  arrayfun(@(n) max(dbfit.(molcalib{i}).all((n),1).r2max),1:nROI)';  
end
w = cat(1,r2fit{:});
w = w/sum(w);
y = cat(1,areacalib{:});
x = cat(1,molH{:});
[b,bint] = regress(w(2:end).*y(2:end),w(2:end).*x(2:end));
dbcalibration.calibcommon.b = b;
dbcalibration.calibcommon.bint = bint;

%% DBCORRSUB : database of multiple correlation between substances (nmrloadmatrix)
dbcorrname = 'dbcorrelation_substances_full.mat';
if exist(fullfile(local,dbcorrname),'file')
    dbcorrsub = load(fullfile(local,'dbcorrelation_substances_full.mat'));
else
    dbcorrsubfull = nmrloadmatrix(fitmol,fitmol,dbpur,dbxpur,dbfit,'verbose',2); % matrix of correlation between substances
    dbcorrsubfullquanti = nmrloadmatrix(fitmol,fitmol,dbpur,dbxpur,dbfit,'verbose',2,'quanti','dbcalibration',dbcalibration,'normvalue',1); % quanti
    save(fullfile(local,dbcorrname),'dbcorrsubfull','dbcorrsubfullquanti');
end

%% DBCLASS : database of classification of susbtance NMR spectra from pairwise correlation coefficient
dbclassname = 'dendrogram_T_worse_full.mat';
if exist(fullfile(local,dbclassname),'file')
   dbclass = load(fullfile(local,dbclassname));
else
    % calcul of pairwise correlation coefficient
    dbconfusionname = 'confusion_matrix_worse_full.mat';  
    if exist(fullfile(local,dbconfusionname),'file')
        load(fullfile(local,dbconfusionname));
    else
        dbconfusion = nmrloadcorrcoeff(mol,[],dbpur,dbxpur,dbfit,'worse',10);
        save(fullfile(local,dbconfusionname),'dbconfusion')
    end
    % calculation , equation
    corrmean = zeros(nmol,nmol);
    corrmeanweight = zeros(nmol,nmol);
    for i= 1:nmol
       for j = 1:nmol
            corrmean(i,j) = dbconfusion.corrmean{i,j}(1,3); 
            corrmeanweight(i,j) = dbconfusion.corrmeanweight{i,j}(1,3); 
       end 
    end
    d2 = @(rho2) 1./(rho2+1e-3)-1; % equation to transform rho2 into distance
    %  ------------ CORRMEAN confusion matrix with normal mean of corr.coeff------------------
    recalcmatrix = sqrt(corrmean.*(corrmean)');
    dstrecalcmatrix = d2(recalcmatrix);
    for i = 1:nmol, for j = i, dstrecalcmatrix(i,j) = 0; end, end
    % cluster data
    X = squareform(dstrecalcmatrix); 
    Z = linkage(X,'average');
    classtotest = 1:nmol;
    T = zeros(nmol,nmol);
    rhobyclass = zeros(nmol,nmol);
    nmolbyclass = zeros(nmol,nmol);
    rhototest = zeros(length(classtotest),5); %rho of each classtotest [min,25prct,mean,75prct,max]
    nmoltotest = zeros(length(classtotest),3); %number mol in each classtotest [µin, mean, max]
    for i = 1:nmol
        T(:,i) = cluster(Z,'maxclust',classtotest(i)); 
        % mol, rho  for each class
        Tunq = unique(T(:,i));
        for iclass = 1: length(Tunq);
            currentclass = find(T(:,i)==Tunq(iclass));
            tmp = recalcmatrix(currentclass,currentclass); % rho2 in current classes 
            rhobyclass(iclass,i) = mean(tmp(:)); % mean of rho2 in current classes 
            nmolbyclass(iclass,i) = length(currentclass); % number of mol in currentclass
        end
        % summary rho fot each classtotest
        rtmp = rhobyclass(:,i);
        nmoltmp = nmolbyclass(:,i);
        rhototest(i,:) = [min(rtmp(rtmp>0)) prctile(rtmp(rtmp>0),25) mean(rtmp(rtmp>0)) prctile(rtmp(rtmp>0),75) max(rtmp(rtmp>0))];
        nmoltotest(i,:) = [min(nmoltmp(nmoltmp>0)) mean(nmoltmp(nmoltmp>0)) max(nmoltmp(nmoltmp>0))];
    end
    % filt & threshold 
    rhofilt = filtzero(rhototest,4); % filt rho
    nmolfilt = filtzero(nmoltotest,3); % filt number of mol
    % matrix rho2 mean by class for all substances -> if 2 substances are in the same class, they have a same rho2mean
    rhomeanclass = zeros(nmol,nmol); 
    for i = 1:nmol
        for j = 1:length(classtotest)
            Tunq = unique(T(:,j));
            for iclass = 1:length(Tunq)
                if T(i,j)==Tunq(iclass), rhomeanclass(i,j) = rhobyclass(iclass,j); end
            end
        end
    end

    % ----corrmeanweight : for confusion matrix with weighted mean of corr.coeff------------------------%
    recalcmatrixw = sqrt(corrmeanweight.*(corrmeanweight)'); %(corrmean3w + (corrmean3w)')./2;
    dstrecalcmatrixw = d2(recalcmatrixw);
    for i = 1:nmol, for j = i, dstrecalcmatrixw(i,j) = 0; end, end
    % cluster data
    Xw = squareform(dstrecalcmatrixw); 
    Zw = linkage(Xw,'average');
    classtotest = 1:nmol;
    Tw = zeros(nmol,nmol);
    rhowbyclass = zeros(nmol,nmol);
    nmolwbyclass = zeros(nmol,nmol);
    rhowtotest = zeros(length(classtotest),5); %rho of each classtotest [min,25prct,mean,75prct,max]
    nmolwtotest = zeros(length(classtotest),3); %number mol in each classtotest [µin, mean, max]
    for i = 1:nmol
        Tw(:,i) = cluster(Zw,'maxclust',classtotest(i)); 
        % mol, rho  for each class
        Twunq = unique(Tw(:,i));
        for iclass = 1: length(Twunq);
            currentclass = find(Tw(:,i)==Twunq(iclass));
            tmp = recalcmatrixw(currentclass,currentclass); % rho2 in current classes 
            rhowbyclass(iclass,i) = mean(tmp(:)); % mean of rho2 in current classes 
            nmolwbyclass(iclass,i) = length(currentclass); % number of mol in currentclass
        end
        % summary rho fot each classtotest
        rtmp = rhowbyclass(:,i);
        nmoltmp = nmolwbyclass(:,i);
        rhowtotest(i,:) = [min(rtmp(rtmp>0)) prctile(rtmp(rtmp>0),25) mean(rtmp(rtmp>0)) prctile(rtmp(rtmp>0),75) max(rtmp(rtmp>0))];
        nmolwtotest(i,:) = [min(nmoltmp(nmoltmp>0)) mean(nmoltmp(nmoltmp>0)) max(nmoltmp(nmoltmp>0))];
    end

    % filt & threshold 
    rhowfilt = filtzero(rhowtotest,4); % filt rho
    nmolwfilt = filtzero(nmolwtotest,3); % filt number of mol
    % matrix rho2 mean by class for all substances -> if 2 substances are in the same class, they have a same rho2mean
    rhowmeanclass = zeros(nmol,nmol); 
    for i = 1:nmol
        for j = 1:length(classtotest)
            Twunq = unique(Tw(:,j));
            for iclass = 1:length(Twunq)
                if Tw(i,j)==Twunq(iclass), rhowmeanclass(i,j) = rhowbyclass(iclass,j); end
            end
        end
    end
    % ---- SAVE -----
    save(fullfile(local,dbclassname),'Z','T','rhomeanclass','rhofilt','nmolfilt','Zw','Tw','rhowmeanclass','rhowfilt','nmolwfilt')
end

%% general nfo
nfo(1).date = date;
nfo(1).local = localname;

%% save all db in .mat
% nmrdb = struct('dbpur',dbpur,'dbxpur',dbxpur,'dbmask',dbmask,'mask',mask,'dbmultiplet',dbmultiplet,'dbmix',dbmix,'dbxmix',dbxmix,...
%                'dbcalibration',dbcalibration,'dbxcalibration',dbxcalibration,'dbfit',dbfit,'dbclass',dbclass,'dbcorrsub',dbcorrsub,'nfo',nfo);
dbname = 'nmrdatabase.mat';
save(fullfile(rootdir(local),dbname),'dbpur','dbxpur','dbmask','mask','dbmultiplet','dbmix','dbxmix','dbcalibration','dbxcalibration','dbfit','dbclass','dbcorrsub','nfo') 

%% plot calibration common
col = cbrewer('qual','Dark2',nmolcalib);
figname = 'calibration_common';
hfig = figure; formatfig(hfig,'figname',figname,'paperposition',[[0.3397    6.4177   20.3046   16.8421]]);
hold on
for i = 1:nmolcalib
    plot(molH{i},areacalib{i},'o','line','none','markersize',10,'markerfacecolor',col(i,:),'markeredgecolor',col(i,:))    
end
xlabel('concentration in proton (mol of proton\cdotL^{-1})'), ylabel('Area')
hl = legend(molcalib,'Location','EastOutside','fontsize',6); set(hl,'box','off')
formatax(gca,'xscale','log','yscale','log')
plot(x(2:end),b*x(2:end),'k','linewidth',1.2)
% print_pdf(300,figname,figurefolder,'nocheck','PaperPositionMode','auto')
% print_png(300,figname,figurefolder,'',0,0,0)

%% plot
figurefolder = fullfile(find_path_toolbox('rmnspec'),'figures','calibration'); if ~exist(figurefolder,'dir'), mkdir(figurefolder); end 

nrows = ceil(sqrt(nmol));
ncols = ceil(nmol/nrows);
figname = 'calibration';
hfig = figure; formatfig(hfig,'figname',figname,'paperposition',[ 0.3387    0.7089   29.0000   19.5663],'paperorientation','landscape');
hs = subplots(ones(1,ncols),ones(1,nrows),0.05,0.05); 
for i = 1:nmol
    subplot(hs(i)), plot(dbcalibration(1).(mol{i}).H,dbcalibration(1).(mol{i}).area,'o')
    hold on, plot(dbcalibration(1).(mol{i}).H,polyval(dbcalibration(1).(mol{i}).p,dbcalibration(1).(mol{i}).H))
    axis tight
    limx = xlim; limy = ylim;
    [irow,icol] = ind2sub([nrows ncols],i);
    if irow==nrows, xlabel('number of proton','fontsize',10), end
    if icol==1, ylabel('Area','fontsize',10), end
    titles(hs(i),mol{i},'fontsize',8)
    formatax(hs(i),'ytick',linspace(min(limy),max(limy),3),'xtick',ceil(linspace(min(limx),max(limx+1),3)),'fontsize',8)
end
set(hs(i+1:end),'visible','off')
% print_pdf(300,figname,figurefolder,'nocheck','PaperPositionMode','auto')
% print_png(300,figname,figurefolder,'',0,0,0)

%% plot
figurefolder = fullfile(find_path_toolbox('rmnspec'),'figures','calibration'); if ~exist(figurefolder,'dir'), mkdir(figurefolder); end 

nrows = ceil(sqrt(nmol));
ncols = ceil(nmol/nrows);
figname = 'calibration_molproton';
hfig = figure; formatfig(hfig,'figname',figname,'paperposition',[ 0.3387    0.7089   29.0000   19.5663],'paperorientation','landscape');
hs = subplots(ones(1,ncols),ones(1,nrows),0.05,0.05); 
for i = 1:nmol
    subplot(hs(i)), plot(dbcalibration(1).(mol{i}).molH,dbcalibration(1).(mol{i}).area*dbpur.(mol{i}).normvalue,'o')
    hold on, plot(dbcalibration(1).(mol{i}).molH,polyval(dbcalibration(1).(mol{i}).pmolproton,dbcalibration(1).(mol{i}).molH))
    axis tight
    limx = xlim; limy = ylim;
    [irow,icol] = ind2sub([nrows ncols],i);
    if irow==nrows, xlabel('mol of proton','fontsize',10), end
    if icol==1, ylabel('Area','fontsize',10), end
    titles(hs(i),mol{i},'fontsize',8)
%     formatax(hs(i),'ytick',linspace(min(limy),max(limy),3),'xtick',ceil(linspace(min(limx),max(limx+1),3)),'fontsize',8)
end
set(hs(i+1:end),'visible','off')
% print_pdf(300,figname,figurefolder,'nocheck','PaperPositionMode','auto')
% print_png(300,figname,figurefolder,'',0,0,0)

%% check k = A/(sum(H)*C(mol/l) for all substances
area = cellfun(@(m) sum(dbcalibration.(m).area),mol);
proton = cellfun(@(m) sum(dbcalibration.(m).H),mol); % from expert, mask file
proton2 = str2double(cellfun(@(m) uncell(regexp(regexp(dbpur.(m).formula,'H\d+','match'),'\d+','match')),mol)); % from formula
conc = cellfun(@(m) (dbpur.(m).concentration*1e-3)/dbpur.(m).Mw,mol); %conversion mg/l to mol/l
norm1 = area./(proton2.*conc); norm1(norm1==Inf) = 0;
norm2 = area./proton2; 
%plot
figname = 'normalisation';
hfig = figure; formatfig(hfig,'figname',figname,'paperposition',[1.5887    2.8778   26.5000   15.2284])
hs = subplots(1,[.9 .1],0,0,'alive',1);
subplot(hs(1)), plot(norm1,'o')
ylabel('A/(H*m)','fontsize',14)
formatax(gca,'yscale','log','xtick',0:1:nmol+1,'xticklabel','','xlim',[0 nmol+1],'fontsize',14)
for i = 1:nmol, text(i,min(ylim)*.9,mol{i},'fontsize',9,'VerticalAlignment','middle','HorizontalAlignment','right','rotation',90), end
% print_pdf(300,figname,figurefolder,'nocheck','PaperPositionMode','auto')
% print_png(300,figname,figurefolder,'',0,0,0)
figname = 'hist_normalisation';
hfig = figure; formatfig(hfig,'figname',figname,'paperposition',[1.5887    2.8778   26.5000   15.2284])
hs = subplots(1,[.9 .1],0,0,'alive',1);
subplot(hs(1)), [~,max1,med1] = histfitlog(norm1);

%plot
figname = 'normalisation';
hfig = figure; formatfig(hfig,'figname',figname,'paperposition',[1.5887    2.8778   26.5000   15.2284])
hs = subplots(1,[.9 .1],0,0,'alive',1);
subplot(hs(1)), plot(norm2,'o')
ylabel('A/H','fontsize',14)
formatax(gca,'yscale','log','xtick',0:1:nmol+1,'xticklabel','','xlim',[0 nmol+1],'fontsize',14)
for i = 1:nmol, text(i,min(ylim)*.9,mol{i},'fontsize',9,'VerticalAlignment','middle','HorizontalAlignment','right','rotation',90), end
% print_pdf(300,figname,figurefolder,'nocheck','PaperPositionMode','auto')
% print_png(300,figname,figurefolder,'',0,0,0)
figname = 'hist_normalisation';
hfig = figure; formatfig(hfig,'figname',figname,'paperposition',[1.5887    2.8778   26.5000   15.2284])
hs = subplots(1,[.9 .1],0,0,'alive',1);
subplot(hs(1)), [~,max2,med2] = histfitlog(norm2);
