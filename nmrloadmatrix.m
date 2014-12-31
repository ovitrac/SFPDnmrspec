function dbout = nmrloadmatrix(listsubstancetest,listsubstanceref,dbpur,dbxpur,dbfit,varargin)
%NMRLOADMATRIX builds matrix of correlation between substances, between substance and mixtures in NMR databases
% syntax: dbout = nmrloadmatrix(listsubstancetest,listsubstanceref,dbxpur,dbfit)
% INPUTS
% listsubstancetest: 1 x m cell {'mol1' 'mol2' 'mol3'}, substances to test (theoretical, fitted spectra), in row
%  listsubstanceref: 1 x n cell {'mola' 'molb'} susbtances building a mixture (if empty = {listsubstancetest}) , in column
%                    or m x n array of mixture spectra 
%             dbpur: nmr database containing real substance spectra data
%            dbxpur: nmr database containing real substance spectra data
%             dbfit: nmr database containing theoretical substance spectra data
% VARARGIN
%       idsubstance: types of keys for link with substances ('commonname' (default), 'reference')
%         dbfitpath: location of fitting NMR spectra database (default = fullfile(find_path_toolbox('rmnspec'),'data_pur'))
%         dbfitname: fitting NMR spectra database (default = dbfit.mat, 
%         normvalue: normalized value of intensity (default = 1e-8)
%           verbose: level for details of dbout [0 1 2 3] (default = 0)
%         proximity: n x 1 array, rang of susbtances by proximity (by calculation of clustering) (default = [], no ranking by proximity)
%            maxlag: array, value maxima of lags (in ppm) to define max.corr. default = [], empty -> consider a whole segment)
% OUTPUTS 
% -------------------- old version (before 02/09/13) ----------------------
%         DBOUT: structure with fields
% VERBOSE = 0 (default)
%           corrmax: k x n matrix containing maxima of correlation 
%                    k: total number of motifs (bands, ROIs) of theoretical database
%                    n: number of substances/ spectra
%       corrzerolag: k x n matrix containing correlation at zero lag
%              rho2: k x n matrix containing correlation coefficients (calculated by NMRCORRCOEFF)
%       subtestlist: m x 1 cell, list of output substances test
%        subreflist: n x 1 cell, list of output substance reference (by alphabetic order)
% VERBOSE = 1
% subtestlistdetail: k x 1 cell, list of substances test with accounting number of motifs (bands, ROIs) for each subtance
%              gate: k x 2 array, [gatemin gatemax] to identify gate of each motif (band, ROI)  
%        lengthtest: k x 1 array, length of each motif test (band, ROI)
%              imax: k x n array, index of maxima correlation
%        centerpeak: k x n array, center position of band in ppm unit
% VERBOSE = 2
%              corr: k x n cell, correlation function  
%              lags: k x n cell, lags when appy correlation function 
%          rho2corr: k x n structures, structs with fields like CORR (3th argument output) of NMRCORRCOEFF 
% VERBOSE = 3
%             Itest: k x n cell, itensity of motifs test (bands, ROIs) (fitting pattern)
%              Iref: k x n cell, itensity of corresponding fragment of reference spectra
%
% ------------------ new structure of DBOUT -------------------------------
%  DBOUT: m x n structure with fields
% VERBOSE = 0 (default)
%           corrmax: k x 1 array containing maxima of correlation 
%                    k: number of motifs (bands, ROIs) of substance
%       corrzerolag: k x 1 array containing correlation at zero lag
%              rho2: k x 1 array containing correlation coefficients (calculated by NMRCORRCOEFF)
%       subtestlist: 1 cell, list of output substances test
%        subreflist: 1 cell, list of output substance reference (by alphabetic order)
%            weight: k x 1 array containing weight attributed by expert in dbmask
% VERBOSE = 1
% subtestlistdetail: 1 cell, list of substances test with accounting number of motifs (bands, ROIs) for each subtance
%              gate: k x 2 array, [gatemin gatemax] to identify gate of each motif (band, ROI)  
%        lengthtest: k x 1 array, length of each motif test (band, ROI)
%              imax: k x 1 array, index of maxima correlation
%        centerpeak: k x 1 array, center position of band in ppm unit
% VERBOSE = 2
%              corr: k x 1 cell, correlation function  
%              lags: k x 1 cell, lags when appy correlation function           
% VERBOSE = 3
%          rho2corr: k x 1 structures, structs with fields like CORR (3th argument output) of NMRCORRCOEFF 
%             Itest: k x 1 cell, itensity of motifs test (bands, ROIs) (fitting pattern)
%              Iref: k x 1 cell, itensity of corresponding fragment of reference spectra
%
% EXAMPLE 1: deconvolution of reference mixture
% dbfit = nmrloaddbfit('path',fullfile(find_path_toolbox('rmnspec'),'data_pur'),'dbname','dbfit.mat');
% [dbpur,dbxpur,dbmask,mask,dbmultiplet,dbmix,dbxmix] = nmrloadbase;
% load(fullfile(find_path_toolbox('rmnspec'),'data_pur','confusion_matrix_worse.mat')); %confusion matrix
% corrmeanweight = zeros(nfitmol,nfitmol);
% for i= 1:nfitmol, for j = 1:nfitmol, corrmeanweight(i,j) = dbconfusion.corrmeanweight{i,j}(1,3); end, end
% confmat = (corrmeanweight + (corrmeanweight)')./2;
% d2 = @(rho2) 1./sqrt(rho2+1e-3)-1; dstmat = d2(confmat);
% for i = 1:length(fieldnames(dbfit)), for j = i, dstmat(i,j) = 0; end, end
% X = squareform(dstmat); Z = linkage(X,'average');
% figure, [~,~,rang] = dendrogram(Z,0); 
% dbcorrsub = nmrloadmatrix({'Oleamide' 'Erucamide' 'Chimassorb944' 'Tinuvin622'},[],dbpur,dbxpur,dbfit,'proximity',rang,'verbose',2,'normvalue',5,'quanti','dbcalibration',dbcalibration);
% dbcorrmix = nmrloadmatrix({'Oleamide' 'Erucamide' 'Chimassorb944' 'Tinuvin622'},[dbmix.N1413.In dbmix.N1412.In],dbxpur,dbfit,'proximity',rang);
% W = lsqlin(dbcorrsub.corrzerolag,dbcorrmix.corrzerolag,[],[],[],[],0);
% 
% EXAMPLE 2: deconvolution of theoretical (artificial) mixture
% dbfit = nmrloaddbfit('path',fullfile(find_path_toolbox('rmnspec'),'data_pur'),'dbname','dbfit.mat');
% [~,dbxpur] = nmrloadbase;
% dbcorrsub = nmrloadmatrix({'Oleamide' 'Erucamide' 'Chimassorb944' 'Tinuvin622'},[],dbxpur,dbfit,'proximity',rang);
% dbcorrmix = nmrloadmatrix({'Oleamide' 'Erucamide' 'Chimassorb944' 'Tinuvin622'},{'Erucamide' 'Chimassorb944'},dbxpur,dbfit,'proximity',rang);
% sample = [.2;.5];
% Cmaxmix = dbcorrmix.corrmax*sample; Czerolagmix = dbcorrmix.corrzerolag*sample;
% W = lsqlin(dbcorrsub.corrmax,Cmaxmix,[],[],[],[],0);
% 
% EXAMPLE 3: full matrix of NMR database
% dbfit = nmrloaddbfit('path',fullfile(find_path_toolbox('rmnspec'),'data_pur'),'dbname','dbfit.mat');
% [dbpur,dbxpur,dbmask,mask,dbmultiplet,dbmix,dbxmix] = nmrloadbase;
% mol = fieldnames(rmfield(dbfit,'help')); 
% data = load(fullfile(find_path_toolbox('rmnspec'),'data_pur','confusion_matrix_worse.mat')); %confusion matrix
% confmat = (data.corrmean + (data.corrmean)')./2;
% d2 = @(rho2) 1./sqrt(rho2+1e-3)-1; dstmat = d2(confmat);
% for i = 1:length(fieldnames(dbfit)), for j = i, dstmat(i,j) = 0; end, end
% X = squareform(dstmat); Z = linkage(X,'average');
% figure, [~,~,rang] = dendrogram(Z,0); 
% dbcorrsub = nmrloadmatrix(mol,[],dbxpur,dbfit,'proximity',rang,'verbose',3);


% RMNSPEC v 0.1 - 28/06/2013 - INRA\Olivier Vitrac, LNE\Mai Nguyen - rev.05/11/14

% History
% 01/07/13 add plot
% 03/07/13 add value of correlation at zero lag -> 2 matrices: for max correlation and for zerolag
% 05/07/13 change significatively inputs 
% 23/07/13 add calculation of rho2 (use nmrcorrcoeff)
% 02/09/13 change structure of dbout, indexing type, group motifs/bands by substance -> dbout = mxn size
%          add maxlag: maxima of lags to define max.corr
% 01/10/13 add dbpur as INPUT and add field weight in dbcorr (weight attributed by expert in dbmask)
% 05/11/13 add keyword 'quanti' to normalize I by calibration data (slope*totoal number of H)
% 06/11/14 check legnht of yfit and yvalid
% defautl
default = struct('idsubstance','commonname',...
                 'dbfitpath',fullfile(find_path_toolbox('rmnspec'),'data_pur'),...
                 'dbfitname','dbfit.mat',...
                 'normvalue',1e-8,...
                 'verbose',0,...
                 'proximity',[],...
                 'maxlag',[],...
                 'dbcalibration',[]);
keyword = ('quanti');
% arg check
o = argcheck(varargin,default,keyword,'nostructexpand');
if nargin<1, error('1 argument1 are required'), end
if nargin<2, listsubstanceref = {}; end
if nargin<3, dbpur = []; end
if nargin<4, dbxpur = []; end
if nargin<5, dbfit = []; end
if isempty(dbpur), dbpur = nmrloadbase; end
if isempty(dbxpur), [~,dbxpur] = nmrloadbase; end %nmrloadascii('pur','calib');
if isempty(dbfit), dbfit = nmrloaddbfit('path',o.dbfitpath,'dbname',o.dbfitname); end
ppm = dbxpur.ppm; commonname = fieldnames(rmfield(dbfit,'help')); step = dbxpur.step; 
if ~isempty(o.maxlag), maxlag = ceil(o.maxlag./step); end

if o.quanti && isempty(o.dbcalibration)
    o.dbcalibration = nmrloaddbspec('dbcalibration');
    nfo = fileinfo(fullfile(find_path_toolbox('rmnspec'),'nmrdatabase.mat'));
    fprintf('DBCALIBRATION is not defined, use dbcalibration in following base:\n%s\ndate and bytes:%s%s',strcat(nfo.path,nfo.filename),nfo.date,nfo.bytes)
end

% check test substances (row-wise)
%   -> set indexes of substances
%   -> set proximity order if required
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
if ~isempty(o.proximity), rang = o.proximity(:); % ranking by proximity of spectra
    indbyproximity = zeros(1,nlistsubstancetest);
    for i =1:nlistsubstancetest, indbyproximity(1,i) = find(rang==indsubtest(i)); end
    indsubtest = rang(sort(indbyproximity,'ascend'));
end                              
subtestname = commonname(indsubtest); % convert indsubtest into commonname, ordered by alphabetic order or proximity order

% check reference substances (column-wise) - reuse test substances if not provided
%   -> set spectra
%   -> set subrefname (with automatic names if not specified)
if ~isempty(listsubstanceref) % --- list of reference substances or spectra
    if iscellstr(listsubstanceref) % list of reference substances
        tmpdb = nmrsubdb(dbxpur,o.idsubstance,listsubstanceref); 
        subrefname = tmpdb.commonname;
        if o.quanti
            spectra = zeros(size(tmpdb.I,1),size(tmpdb.I,2));
            for i = 1:size(tmpdb.I,2)
                spectra(:,i) = tmpdb.I(:,i)/(o.dbcalibration.dbcalibration.(subrefname{i}).p(1)*sum(o.dbcalibration.dbcalibration.(subrefname{i}).H));  % spectra of n mol ref
            end
        else spectra = tmpdb.I;
        end
        nlistsubstanceref = length(listsubstanceref);
    else                    % spectra are supplied
        spectra = listsubstanceref; 
        nlistsubstanceref = size(spectra,2); 
        subrefname = cell(1,nlistsubstanceref);
        for i=1:nlistsubstanceref, subrefname{i} = sprintf('spectra%02d',i); end % default names
    end
else % --- no list of reference substances
    listsubstanceref = listsubstancetest;
    spectra = [];
    subrefname = subtestname;
    nlistsubstanceref = length(listsubstanceref);
end


% main (filling of matrix by row)
narea = 0; for i=1:nlistsubstancetest, narea = narea + size(dbfit.(subtestname{i}).all,1); end
corrmax = zeros(narea,nlistsubstanceref); imax = zeros(narea,nlistsubstanceref);
corrzerolag = zeros(narea,nlistsubstanceref);
rho2 = cell(narea,nlistsubstanceref);
corr = cell(narea,nlistsubstanceref);
lags = cell(narea,nlistsubstanceref);
subtestnametdetail = cell(narea,1);
npatternbysubtest = zeros(nlistsubstancetest,1);
gate = zeros(narea,2);
lengthtest = zeros(narea,1);
centerpeak = zeros(narea,1);
weight = zeros(narea,1);
yfit = cell(narea,1);
yvalid = cell(narea,nlistsubstanceref);
rho2corr = repmat(struct('yref',[],'yrefnorm',[],'y',[],'ynorm',[],'lags',[],...
                         'xgatemin',[],'xgatemax',[],'rho2',[],'npad',[],'typepad',[],'m',[]),narea,1);
k=0;
for i = 1:nlistsubstancetest; % pattern mol, fitting mol (in row)
    npatternbysubtest(i,1) = size(dbfit.(subtestname{i}).all,1);
    for iROI = 1:npatternbysubtest(i,1)
        k=k+1;
        xfit = dbfit.(subtestname{i}).all(iROI,1).ppm; % ppm of peak
        if o.quanti, yfit{k,1} = dbfit.(subtestname{i}).all(iROI,1).Ifit/(o.dbcalibration.dbcalibration.(subtestname{i}).p(1)*sum(o.dbcalibration.dbcalibration.(subtestname{i}).H));
        else yfit{k,1} = dbfit.(subtestname{i}).all(iROI,1).Ifit;
        end          
        valid = ((ppm>=min(xfit))&(ppm<=max(xfit)));
        % check lengh
        if length(find(valid)) > length(xfit)
            ndiff = length(find(valid)) - length(xfit);
            valid(find(valid,1,'first'):find(valid,1,'first')+ndiff-1) = 0;
        elseif length(find(valid)) < length(xfit)
            ndiff = length(xfit) -length(find(valid));
            valid(find(valid,1,'first')-ndiff:find(valid,1,'first')) = 1;
        end
        subtestnametdetail{k,1} = subtestname{i};
        gate(k,:) = [min(xfit) max(xfit)];
        lengthtest(k,1) = length(yfit{k,1}); % lengt of corr = (lengthtest*2)+1
        centerpeak(k,1) = dbfit.(subtestname{i}).all(iROI,1).window.center; % center of motif in ppm
        weight(k,1) = dbpur.(subtestname{i}).gates(iROI,4);% weight defined by user in dbmask   
        if ~isempty(o.maxlag), 
            if maxlag > lengthtest(k,1)
                maxlagcorr = min(lengthtest(k,1),maxlag); maxlag = maxlagcorr;
            else maxlagcorr = max(lengthtest(k,1),maxlag);
            end
        else maxlagcorr = lengthtest(k,1); maxlag = maxlagcorr;
        end
        for j = 1:nlistsubstanceref %reference substances, in column
            if isempty(spectra), 
                if o.quanti, y = dbxpur.I(:,indsubtest(j))/(o.dbcalibration.dbcalibration.(subtestname{j}).p(1)*sum(o.dbcalibration.dbcalibration.(subtestname{j}).H));
                else y = dbxpur.I(:,indsubtest(j)); 
                end
            else y = spectra(:,j); 
            end
            yvalid{k,j} = y(valid);
            [corr{k,j},lags{k,j}] = xcorr(yfit{k,1}/o.normvalue,yvalid{k,j}/o.normvalue,maxlagcorr,'biased');
            [corrmax(k,j),itmp] = max(corr{k,j}((maxlagcorr-maxlag+1):(maxlagcorr-maxlag+1)+2*maxlag));
            imax(k,j) = maxlagcorr-maxlag+itmp; % conversion itmp in imax (in whole segment)
            zerolag = length(yvalid{k,j})+1;
            corrzerolag(k,j) = corr{k,j}(zerolag,1);
            rholags = nmrfindlags(yfit{k,1},yvalid{k,j});
            [tmp,~,tmpcorr] = nmrcorrcoeff(yfit{k,1},yvalid{k,j},rholags); 
            if length(rholags)>1, [~,indrho] = max(tmp(:,3));
                rho2{k,j} = tmp(indrho,:);
                rho2corr(k,j) = tmpcorr(indrho); 
            else
                rho2{k,j} = tmp;
                rho2corr(k,j) = tmpcorr; 
            end
        end
    end
end

% group data by 'substance'
corrmaxgr = cell(nlistsubstancetest,nlistsubstanceref); % max corr of dbcorrsub grouped by substances
corrzerolaggr = cell(nlistsubstancetest,nlistsubstanceref);  % corr at zerolag of dbcorrsub grouped by substances
rho2gr = cell(nlistsubstancetest,nlistsubstanceref);
subtestnametdetailgr = cell(nlistsubstancetest,1);
gategr = cell(nlistsubstancetest,1);
lengthtestgr = cell(nlistsubstancetest,1);
imaxgr = cell(nlistsubstancetest,nlistsubstanceref); % index of max corr of dbcorrsub grouped by substances
centerpeakgr = cell(nlistsubstancetest,1); % center of band of dbcorrsub grouped by substances
weightgr = cell(nlistsubstancetest,1); % center of band of dbcorrsub grouped by substances
corrgr = cell(nlistsubstancetest,nlistsubstanceref); % correlation of dbcorrsub grouped by substances
laggr = cell(nlistsubstancetest,nlistsubstanceref); % lags of dbcorrsub grouped by substances
rho2corrgr = cell(nlistsubstancetest,nlistsubstanceref); 
yfitgr = cell(nlistsubstancetest,1);
yvalidgr = cell(nlistsubstancetest,nlistsubstanceref);

k = 1;
for i = 1:nlistsubstancetest
    nROI = npatternbysubtest(i); % size(dbfit.(subdbcorrsub.subtestlist{i}).all,1);
    for j = 1:nlistsubstanceref
        corrmaxgr{i,j} = corrmax(k:k+nROI-1,j);
        corrzerolaggr{i,j} = corrzerolag(k:k+nROI-1,j);
        rho2gr{i,j} = rho2(k:k+nROI-1,j);
        subtestnametdetailgr{i,1} = subtestnametdetail{k:k+nROI-1,1};
        gategr{i,1} = gate(k:k+nROI-1,:);
        lengthtestgr{i,1} = lengthtest(k:k+nROI-1,1);
        imaxgr{i,j} = imax(k:k+nROI-1,j);
        centerpeakgr{i,1} = centerpeak(k:k+nROI-1,1);
        weightgr{i,1} = weight(k:k+nROI-1,1);
        corrgr{i,j} = corr(k:k+nROI-1,j);
        laggr{i,j} = lags(k:k+nROI-1,j);
        rho2corrgr{i,j} = rho2corr(k:k+nROI-1,j); 
        yfitgr{i,1} = yfit(k:k+nROI-1,1);
        yvalidgr{i,j} = yvalid(k:k+nROI-1,j);

    end
    k = k+nROI;
end

% outputs
if o.verbose > 0 % verbose 1
    dbout = repmat(struct('corrmax',[],'corrzerolag',[],'rho2',[],'subreflist',[],'subtestlist',[],'weight',[],...
                          'subtestlistdetail',[],'npatternbysubtest',[],'gate',[],'lengthtest',[],'imax',[],'centerpeak',[]),nlistsubstancetest,nlistsubstanceref);
elseif o.verbose > 1 % verbose 2
    dbout = repmat(struct('corrmax',[],'corrzerolag',[],'rho2',[],'subreflist',[],'subtestlist',[],'weight',[],... 
                          'subtestlistdetail',[],'npatternbysubtest',[],'gate',[],'lengthtest',[],'imax',[],'centerpeak',[],...
                          'corr',[],'lags',[]),nlistsubstancetest,nlistsubstanceref); 
elseif o.verbose > 2 % verbose 3
    dbout = repmat(struct('corrmax',[],'corrzerolag',[],'rho2',[],'subreflist',[],'subtestlist',[],'weight',[],... verbose = 0
                          'subtestlistdetail',[],'npatternbysubtest',[],'gate',[],'lengthtest',[],'imax',[],'centerpeak',[],...
                          'corr',[],'lags',[],'rho2corr',[],'yfit',[],'yvalid',[]),nlistsubstancetest,nlistsubstanceref); %verbose 2 and 3
else dbout = repmat(struct('corrmax',[],'corrzerolag',[],'rho2',[],'subreflist',[],'subtestlist',[],'weight',[]),nlistsubstancetest,nlistsubstanceref); %verbose = 0
end
for i = 1:nlistsubstancetest
    for j = 1:nlistsubstanceref
        dbout(i,j).corrmax = corrmaxgr{i,j};
        dbout(i,j).corrzerolag = corrzerolaggr{i,j};
        dbout(i,j).rho2 = rho2gr{i,j};
        dbout(i,j).subreflist = subrefname{j};
        dbout(i,j).subtestlist = subtestname{i,1};
        dbout(i,j).weight = weightgr{i,1};
        if o.verbose > 0  % verbose = 1
            dbout(i,j).subtestlistdetail = subtestnametdetailgr{i,1};
            dbout(i,j).npatternbysubtest = npatternbysubtest(i,1);
            dbout(i,j).gate = gategr{i,1};
            dbout(i,j).lengthtest = lengthtestgr{i,1};
            dbout(i,j).imax = imaxgr{i,j};
            dbout(i,j).centerpeak = centerpeakgr{i,1};
        end
        if o.verbose > 1  % verbose = 2
            dbout(i,j).corr = corrgr{i,j};
            dbout(i,j).lags = laggr{i,j};
        end
        if o.verbose > 2  % verbose = 3
            dbout(i,j).rho2corr = rho2corrgr{i,j};
            dbout(i,j).Itest = yfitgr{i,1};
            dbout(i,j).Iref = yvalidgr{i,j};
        end
    end
end

% dbout = struct('corrmax',[],'corrzerolag',[],'rho2',[],'subtestlist',[],'subreflist',[]); % verbose = 0
% dbout.corrmax = corrmax;
% dbout.corrzerolag = corrzerolag;
% dbout.rho2 = rho2;
% dbout.subtestlist = subtestname;
% dbout.subreflist = subrefname;
% if o.verbose > 0  % verbose = 1
%     dbout.subtestlistdetail = subtestnametdetail;
%     dbout.npatternbysubtest = npatternbysubtest;
%     dbout.gate = gate;
%     dbout.lengthtest = lengthtest;
%     dbout.imax = imax;
%     dbout.centerpeak = centerpeak;
% end
% if o.verbose > 1  % verbose = 2
%     dbout.corr = corr;
%     dbout.lags = lags;
% end
% if o.verbose > 2  % verbose = 3
%     dbout.rho2corr = rho2corr;
%     dbout.Itest = yfit;
%     dbout.Iref = yvalid;
% end

%% plots if no output
if ~nargout
    hfig = figure; formatfig(hfig,'figname','correlation','paperposition',[2.2420    0.1774   16.5000   29.0000])
    nrows = narea; ncols = nlistsubstanceref; 
    hs = subplots(ones(1,ncols),ones(1,nrows),0,0);
    k=0;
    for i = 1:nlistsubstanceref;
        for j = 1:narea            
            k=k+1;
            subplot(hs(k)), hold on
            plot(lags{j,i},corr{j,i})
            plot(lags{j,i}(imax(j,i)),corrmax(j,i),'marker','s','markersize',5,'linestyle','none','markerfacecolor','r','markeredgecolor','r') 
            plot(0,corrzerolag(j,i),'marker','o','markersize',5,'linestyle','none','markerfacecolor','k','markeredgecolor','k')
            [irow,icol] = ind2sub([nrows ncols],k);
            if irow==1, title(subrefname{i}), end
            if irow==nrows, xlabel('lags'), end
            if icol==1 && irow==ceil(nrows/2), ylabel('Correlation'), end   
            if icol>1, set(gca,'yticklabel',' '), end
            if irow<nrows, set(gca,'xticklabel',' '), end
        end
    end
    set(hs,'fontsize',7)
end
