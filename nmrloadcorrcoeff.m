function dbout = nmrloadcorrcoeff(listsubstancetest,listsubstanceref,dbpur,dbxpur,dbfit,varargin)
%NMRLOADCORRCOEFF builds matrix of correlation coefficients between substances, between substance and mixtures in NMR databases
% syntax: dbout = nmrloadcorrcoeff(listsubstancetest,listsubstanceref,dbpur,dbxpur,dbfit)
% INPUTS
% listsubstancetest: 1 x m cell {'mol1' 'mol2' 'mol3'}, substances to test(theoretical, fitted spectra), in row
%  listsubstanceref: 1 x n cell {'mola' 'molb'} susbtances to be compared, in column (if empty = {listsubstancetest})
%                    or m x n array of mixture spectra (
%             dbpur: nmr database containing real substance spectra data
%            dbxpur: concatenated nmr database containing real substance spectra data
%             dbfit: nmr database containing theoretical substance spectra data
% VARARGIN
%       idsubstance: types of keys for link with substances ('commonname' (default), 'reference')
%         dbfitpath: location of fitting NMR spectra database (default = fullfile(find_path_toolbox('rmnspec'),'data_pur'))
%         dbfitname: fitting NMR spectra database (default = dbfit.mat, 
%         threshold: threshold of corr.coeff to make gap when corr.coeff < threshold(default=0.6)
%             worse: 'worse' parameter when calculating corr.coeff. (default=1) (set worse = 10 or more when comparing substances and mixtures) 
% OUTPUTS 
%         DBOUT: structure with fields
%          corrmean: mean of corr.coeff (3 types of rho2: normal, small gate and gate inf)
% 	 corrmeanweight: weighted mean of corr.coeff
%       subtestlist: name of tested substances (fitting, pattern)
%        subreflist: name of ref substances (mixture spectra or substances)
%
% EXAMPLE 1: confusion matrix
% dbfit = nmrloaddbfit('path',fullfile(find_path_toolbox('rmnspec'),'data_pur'),'dbname','dbfit.mat');
% [dbpur,dbxpur] = nmrloadbase;
% dbcorr = nmrloadcorrcoeff(fieldnames(dbfit),[],dbpur,dbxpur,dbfit);
% corrmean = zeros(length(fieldnames(dbfit)),length(fieldnames(dbfit)));
% for i= 1:length(fieldnames(dbfit)), for j = 1:length(fieldnames(dbfit)), corrmean(i,j) = dbcorr.corrmean{i,j}(1,3); end, end
% figure, imagesc(corrmean), colormap(cbrewer('seq','OrRd',length(fieldnames(dbfit))))
%
% EXAMPLE 2: deconvolution of reference mixture
% dbfit = nmrloaddbfit('path',fullfile(find_path_toolbox('rmnspec'),'data_pur'),'dbname','dbfit.mat');
% [dbpur,dbxpur,dbmask,mask,dbmultiplet,dbmix,dbxmix] = nmrloadbase;
% dbcorr = nmrloadcorrcoeff(fieldnames(dbfit),dbmix.SFPDPP3.In,dbpur,dbxpur,dbfit,'worse',10);
%
% RMNSPEC v 0.1 - 28/06/2013 - INRA\Olivier Vitrac, LNE\Mai Nguyen - rev.26/11/14
% history
% 26/11/14: add time counter
% default
default = struct('idsubstance','commonname',...
                 'dbfitpath',fullfile(find_path_toolbox('rmnspec'),'data_pur'),...
                 'dbfitname','dbfit.mat',...
                 'threshold',0.6,...
                 'worse',1,...
                 'weight','weight');

% arg check
o = argcheck(varargin,default);
if nargin<1, error('1 argument is required'), end
if nargin<2, listsubstanceref = {}; end
if nargin<3, dbpur = []; end
if nargin<4, dbxpur = []; end
if nargin<5, dbfit = []; end
if isempty(dbpur), dbpur = nmrloadbase; end
if isempty(dbxpur), [~,dbxpur] = nmrloadbase; end
if isempty(dbfit), dbfit = nmrloaddbfit('path',o.dbfitpath,'dbname',o.dbfitname); end
ppm = dbxpur.ppm; commonname = fieldnames(dbfit);

% check test substances (row-wise)
%   -> set indexes of substances
if ~isfield(dbxpur,o.idsubstance), error('the field ''%s'' does not exist in the supplied database',o.idsubstance), end
if ~ischar(dbxpur.(o.idsubstance)) && ~iscellstr(dbxpur.(o.idsubstance)), error('the field ''%s'' in the supplied database must contains only strings',o.idsubstance), end
if ~iscell(listsubstancetest), listsubstancetest = {listsubstancetest}; end
if ~iscellstr(listsubstancetest), error('the substance keys must be given in a cell array of strings (numerical values are excluded)'), end 

nlistsubstancetest = length(listsubstancetest); 
if nlistsubstancetest>length(unique(listsubstancetest)), error('substance keys must be unique, please check'), end
[~,indsubtest] = intersect(dbxpur.(o.idsubstance),listsubstancetest); % indices of test substances in alphabetic order =/ than order in dbxpur
if length(indsubtest)<nlistsubstancetest
    dispf('ERROR\t%d values given as substance keys cannot be found in ''%s''',nlistsubstancetest-length(indsubtest),o.idsubstance)
    cellfun(@(m) dispf('\t''%s'' is missing',m),setdiff(listsubstancetest,dbxpur.(o.idsubstance)))
    error('%d missing values (see above)',nlistsubstancetest-length(indsubtest))
end                            
subtestname = commonname(indsubtest); % convert indsubtest into commonname, ordered by alphabetic order or proximity order

% check reference substances (column-wise) - reuse test substances if not provided
%   -> set spectra
%   -> set subrefname (with automatic names if not specified)
if ~isempty(listsubstanceref) % --- list of reference substances or spectra
    if iscellstr(listsubstanceref) % list of reference substances
        tmpdb = nmrsubdb(dbxpur,o.idsubstance,listsubstanceref); 
        spectra = tmpdb.I;         % spectra of n mol ref
        subrefname = tmpdb.commonname;
        nlistsubstanceref = length(listsubstanceref);
    else                           % spectra are supplied
        spectra = listsubstanceref; 
        nlistsubstanceref = size(spectra,2); 
        subrefname = cell(1,nlistsubstanceref);
        for i=1:nlistsubstanceref, subrefname{i} = sprintf('spectra%02d',i); end % default names
    end
else                          % --- no list of reference substances
    listsubstanceref = listsubstancetest;
    spectra = [];
    subrefname = subtestname;
    nlistsubstanceref = length(listsubstanceref);
end

% main 
corrmean = cell(nlistsubstancetest,nlistsubstanceref);
corrmeanweight = cell(nlistsubstancetest,nlistsubstanceref);
iter = struct('current',0,'total',nlistsubstanceref*nlistsubstancetest,'t0',clock,'executed',0); screen = '';
for i = 1:nlistsubstanceref % in column
    if isempty(spectra), y = dbxpur.I(:,indsubtest(i)); else y = spectra(:,i); end           
    for j = 1:nlistsubstancetest % in row, fitting, pattern
        [ccum,ccumw] = deal(0);
        nROI = size(dbfit.(subtestname{j}).all,1);
        ROI = dbpur.(subtestname{j}).gates; ROI = ROI(ROI(:,4)>0,:);
        weight = ROI(:,4); % weight defined by user in dbmask to calculate weighted mean of corr.coeff.
%         weight = dbpur.(subtestname{j}).(o.weight)(1:nROI);
%         if iscell(weight), weight = str2double(uncell(regexp(weight,'\d+','match'))); end
        w = ones(nROI,3);
        for iROI = 1:nROI;
            xfit = dbfit.(subtestname{j}).all(iROI,1).ppm; %ppm of peak
            yfit = dbfit.(subtestname{j}).all(iROI,1).Ifit;  
            valid = ((ppm>=min(xfit))&(ppm<=max(xfit)));
            yvalid = y(valid);
            lags = nmrfindlags(yfit,yvalid);
            tmp = nmrcorrcoeff(yfit,yvalid,lags); 
            if length(lags)>1, tmp = max(tmp); end            
            w(iROI,tmp<o.threshold) = o.worse;
            ccum = ccum + tmp.*w(iROI,:); % 3 types of rho2
            ccumw = ccumw + tmp.*(w(iROI,:)*weight(iROI)); % 3 types of rho2 + weight
        end
        % display
        iter.current = iter.current + 1;
        iter.executed = iter.current/iter.total;
        dt = etime(clock,iter.t0);
        screen = dispb(screen,'Calculation procedure %0.3g%% - elapsed time %0.4gs - remaining time %0.4gs',iter.executed*100,dt,(dt/iter.executed)*(1-iter.executed));
        % end display
        corrmean{j,i} = ccum./sum(w,1);
        corrmeanweight{j,i} = ccumw./(w'*weight)';
    end
end

% outputs
dbout = struct('corrmean',[],'corrmeanweight',[],'subtestlist',[],'subreflist',[]);
dbout.corrmean = corrmean;
dbout.corrmeanweight = corrmeanweight;
dbout.subtestlist = subtestname;
dbout.subreflist = subrefname;

%% plot
if ~nargout
    corrmean3 = zeros(nlistsubstancetest,nlistsubstanceref);
    corrmean3w = zeros(nlistsubstancetest,nlistsubstanceref);
    for i= 1:nlistsubstancetest 
        for j = 1:nlistsubstanceref
            corrmean3(i,j) = corrmean{i,j}(1,3); corrmean3w(i,j) = corrmeanweight{i,j}(1,3); % rho 2 gate inf
        end
    end
    col = rgb({'Blue' 'Crimson'});
    for i = 1:nlistsubstanceref
        [corrmeansort,isort] = sort(corrmean3(:,i),'descend');
        [corrmeanweightsort,iwsort] = sort(corrmean3w(:,i),'descend'); 
        figname = ['corr' subrefname{i}];
        hfig = figure; formatfig(hfig,'figname',figname,'paperposition',[0.3387    0.7975   29.0000   19.3890],'paperorientation','landscape');
        hs = subplots([.15 .7 .15],[.1 .9],0,0,'alive',4);
        subplot(hs(1)), plot(corrmeansort,1:nlistsubstancetest,'-o','color',col(1,:))
        xlabel('\rho^2_{gate}','fontsize',14)
        ax1 = gca;
        formatax(ax1,'box','off','ylim',[0 nlistsubstancetest+1],'ytick',0:nlistsubstancetest+1,'xlim',[0 1],'xtick',0:0.1:1,'Xcolor',col(1,:),'Ycolor',col(1,:),'yticklabel','')
        ax2 = axes('position',get(ax1,'position'),'xAxisLocation','top','yAxisLocation','right','color','none','XColor',col(2,:),'YColor',col(2,:)); 
        hold on
        plot(corrmeanweightsort,1:nlistsubstancetest,'o-','parent',ax2,'linewidth',1,'color',col(2,:)) 
        formatax(ax2,'box','off','ylim',[0 nlistsubstancetest+1],'ytick',0:nlistsubstancetest+1,'xlim',[0 1],'xtick',0:0.1:1,'yticklabel','') 
        xlabel('\rho^2_{gate,weighted}','fontsize',14)
        title(subrefname{i},'fontsize',14)
        for j=1:nlistsubstancetest
            text(0,j,{'\color{blue}' subtestname{isort(j)}},'fontsize',9,'VerticalAlignment','middle','HorizontalAlignment','right')
            text(1,j,{'\color{red}' subtestname{iwsort(j)}},'fontsize',9,'VerticalAlignment','middle','HorizontalAlignment','left')
        end
    end
end