function R = deformulateNMR(nmrdb,mixturespec,varargin)
%DEFORMUALTENMR deconvolves a 1H NMR spectrum according to a dictionary of 1H NMR spectra using the method [1]
%   R = deformulateNMR(nmrdb,mixture,'property1',value1,'property2',value2,...)
%   listofproperties = deformulateNMR('default') retrieves optional properties
%   deformulateNMR help or deformulateNMR('help') displays help
%
% INPUTS
%              nmrdb: structure created with script "calibration.m" or filename where the structure can be loaded
%            mixture: n x 2 array of mixture spectrum (1st column: chemical shift and 2nd column: intensity)
%                     or txt  file containing 2 columns and without header    
%
% Optional pair property/value (default value)
%            'worse': numeric, 'worse' parameter when calculating corr.coeff. (default=10)
%      'idsubstance': string, type of subtance names in database for calulating corr. coeff (default='commonname')
%        'threshold': numeric, threshol of corr.coeff to make gap when corr.coeff < threshold(default=0.6) 
%       'outputpath': string, directory of output for saving results database 
%      'mixturename': string, name of mixture, for legend (default = 'mixtureS1')
%       'nmolchoose': numeric, number of substances choosen for deconvolution (default=10)
%  'nclasstheorical': numeric, number of classes in spectral database of substances (default=[])
%'rhothresholdclass': numeric, threshold of corr.coeff to sort substances by class (default=0.75)
%         'nmoltree': number substances to be plotted in scenario tree (default=5)
%         'nplotfit': numeric, number of substances to be fitted and ploted
%    'colormapgraph': colormap for tree graph (default=jet(64) in log scale)
%
% OUTPUT
%   R = structure with fields
%
% REFERENCE
% [1]
%
% PRINCIPLE OF THE METHOD
% List of tasks carried out in DEFORMUALTENMR
% step 1: pair correlation coefficient between spectral database and mixture (nmrloadcorrcoeff)
% step 2: extract the best "10" substances or rho > 0.75 or ... 
% step 3: multiple correlation (nmrloadmatrix)
% step 4: inversion of matrix (nmrlsqnonngeg)
%
% See also: nmrlsqnonneg, nmrloadbase, nmrloadascii, nmrloaddbfit, nmrloadmatrix, nmrloadcorrcoeff, plotnmrlsqnonneg

% RMNSPEC v 0.5 - 30/09/2014 - INRA\Olivier Vitrac, LNE\Mai Nguyen - rev. 15/10/2014
% History
% 15/10/2014: use lacalname only if outputpathdefault is empty
%
% default
outputpathdefault = '';
default = struct('worse',10,...
                 'idsubstance','commonname',... 
                 'threshold',0.6,... 
                 'outputpath',outputpathdefault,...
                 'mixturename','mixtureS1',...
                 'nmolchoose',10,...
                 'nclasstheorical',[],...
                 'rhothresholdclass',0.75,...
                 'nmoltree',5,...
                 'nplotfit',4,...
                 'colormapgraph',[]);

% help
main = struct('varname',{'nmrdb';'mixture'},'help',...
                  {'structure created with script "calibration.m" or filename where the structure can be loaded'
                   'n x 2 array of mixture spectrum (1st column: chemical shift and 2nd column: intensity) or txt file containing 2 columns and without header' });
default.help.worse = 'numeric, worse parameter when calculating correlation coefficient (default=10)';
default.help.idsubstance = 'string, type of subtance names in database for calulating correlation coefficient (default=commonname)';
default.help.threshold = 'numeric, threshol of correlation coefficient to make gap when corr.coeff < threshold (default=0.6)';
default.help.outputpath = 'string, directory of output for saving results database';
default.help.mixturename = 'string, name of mixture, for legend (default = mixtureS1)';
default.help.nmolchoose = 'numeric, number of substances choosen for deconvolution (default=10)';
default.help.nmoltree = 'number substances to be plotted in scenario tree (default=5)';
default.help.nplotfit = 'numeric, number of substances to be fitted and ploted';
default.help.colormapgraph ='colormap for tree graph (default=jet(64)';
default.help.nclasstheorical = 'numeric, number of classes in spectral database of substances (default=[])';
default.help.rhothresholdclass = 'numeric, threshold of corr.coeff to sort substances by class (default=0.75)';

% COMMAND SYNTAXES (help, default)
if (nargin==1) && ischar(varargin{1})
    suffix = {'st' 'nd' 'rd' 'th'};
    if strcmpi(varargin{1},'help')
        % help syntax: deformulateNMR help or deformulateNMR('help')
        dispf('\n\n\tdeformulateNMR help\n')
        dispf('\tR = deformulateNMR(nmrdb,mixture,''property1'',value1,''property2'',value2,...)\n')
        dispf('\t list of mandatory pair property/value (default value)')
        for i=1:length(main), dispf('[%d%s argument] %10s: %s',i,suffix{i},main(i).varname,main(i).help), end
        dispf('\n\t list of optional pair property/value (default value)')
        for prop = sort(fieldnames(default.help))'; dispf('%20s: %s',prop{1},default.help.(prop{1})), end
        if nargout, R = default.help; end
    elseif strcmpi(varargin{1},'default')
        % default syntax: R = deformulateNMR('default')
        R = rmfield(default,'help');
    end
    return
end

% ARGCHECK
o = argcheck(varargin,default,'nostructexpand','case');
% MANDATORY PROPERTIES: nmrdb,mixture
if nargin < 2, error('2 arguments are required'), end
dbinfo = '';
if ischar(nmrdb)
    dbinfo = fileinfo(nmrdb);
    if exist(nmrdb,'file')
        dispf('\tload reference NMR database:'), fileinfo(nmrdb), load(nmrdb)
    else
        error('the reference NMR database ''%s'' does not exist in ''%s''',lastdir(nmrdb),rootdir(nmrdb))
    end
end
if ~isstruct(nmrdb) || ~isfield(nmrdb,'dbpur') || ~isfield(nmrdb,'dbxpur') || ~isfield(nmrdb,'dbfit') || ~isfield(nmrdb,'dbclass') || ~isfield(nmrdb,'dbcalibration') || ~isfield(nmrdb,'dbcorrsub')
    error('nmrdb must be created with script ''calibration.m''')
end
if ischar(mixturespec)
    if exist(mixturespec,'file')
        dispf('\tload reference NMR database:'), fileinfo(mixturespec), mixturespec = load(mixturespec);
    else
        error('the mixture spectum ''%s'' does not exist in ''%s''',lastdir(mixturespec),rootdir(mixturespec))
    end
end
if ~isnumeric(mixturespec) || size(mixturespec,2)~=2, error('mixture must be nx2 numeric vector (1st column: ppm and 2nd column: I'), end
if ~isnumeric(mixturespec) || size(mixturespec,1)~=length(nmrdb.dbxpur.ppm) || mixturespec(1,1)~=nmrdb.dbxpur.ppm(1) || mixturespec(end,1)~=nmrdb.dbxpur.ppm(end)
    error(sprintf('Please define mixture NMR spectrum with chemical shift limited in %d-%d ppm and containing %d points',nmrdb.dbxpur.ppm(1),nmrdb.dbxpur.ppm(end),length(nmrdb.dbxpur.ppm)))
end
if isempty(o.outputpath)
    switch localname
        case 'LP-MOL5'
            o.outputpath = 'C:\Data\Olivier\INRA\Projects\SafeFoodPack_design\livrables\D1_deformulation_rule';
        case 'LP_MOL2'
            o.outputpath = 'C:\Data\Olivier\INRA\Projects\SafeFoodPack_design\livrables\D1_deformulation_rule';      
    end
end
% extract substance names for calculation : commonname
mol = fieldnames(rmfield(nmrdb.dbfit,'help'));
repet = cellfun(@isempty,regexp(mol,'R\d+'));% extract only fitmol without repetition
mol = mol(repet); nmol = length(mol);
% colormap
if ~isempty(o.colormapgraph), colormapgraph = o.colormapgraph; 
else % log color sclae
    ncolor = 64; xcolor = linspace(0,1,ncolor);
    cmap = jet(ncolor);
    funccolor = @(x,a) log(1+a*x);
    cmaplog = @(a) interp1(xcolor,cmap,(funccolor(xcolor,a)-funccolor(xcolor(1),a))/(funccolor(xcolor(end),a)-funccolor(xcolor(1),a)));
    colormapgraph = cmaplog(50);
end
% classification, number of classes
if ~isempty(o.nclasstheorical), nclasstheorical = o.nclasstheorical;
elseif ~isempty(o.rhothresholdclass), 
    nclasstheorical = ceil(csaps(nmrdb.dbclass.rhowfilt(:,3),1:size(nmrdb.dbclass.rhowfilt,1),[],o.rhothresholdclass)); %rhofilt(:,3) = mean of rho2 %number of classes to be tested
else warning(sprintf('No information about the classification of spectral database. Number of classes il the number of substances: %d',nmol))
end
% output directory
outputpath = fullfile(o.outputpath,o.mixturename);
if ~exist(outputpath,'dir'), mkdir(outputpath); end 

%% STEP 1: pairwise correlation coefficient 
dbcorrcoeffname = sprintf('dbcorrcoeff_mixture_%s.mat',o.mixturename);
if exist(fullfile(outputpath,dbcorrcoeffname),'file')
    load(fullfile(outputpath,dbcorrcoeffname));
else dbcorrcoeff = nmrloadcorrcoeff(mol,mixturespec(:,2),nmrdb.dbpur,nmrdb.dbxpur,nmrdb.dbfit,'worse',o.worse,'idsubstance',o.idsubstance,'threshold',o.threshold);      
    if ~exist(outputpath,'dir'), mkdir(outputpath); end 
    save(fullfile(outputpath,dbcorrcoeffname),'dbcorrcoeff')
end     

%% STEP 2: choose of substances for deconvolution 
corrmeanw = cell2mat(arrayfun(@(m) dbcorrcoeff.corrmeanweight{(m),1}(1,3),1:nmol,'UniformOutput',false))';
[corrmeanweightsort,iwsort] = sort(corrmeanw,'descend'); %#ok<UDIM>
molsortbyrho = dbcorrcoeff.subtestlist(iwsort);

tmpmol = molsortbyrho(1:o.nmolchoose);
tmpmol = setdiff(tmpmol,{'PE' 'PP' 'PS'}');
lengthtest = o.nmolchoose+o.nmolchoose-length(tmpmol);
while length(tmpmol) < o.nmolchoose
    tmpmol = molsortbyrho(1:lengthtest);
    tmpmol = setdiff(tmpmol,{'PE' 'PP' 'PS'}');
    lengthtest = lengthtest+1;
end    
moltestsortbyrho = molsortbyrho(ismember(molsortbyrho,tmpmol)); % list of n mols test in decreasing order of rho2
rhotestsort = corrmeanweightsort((ismember(molsortbyrho,tmpmol)),1);% list of rho2 of n mols test
 
%% STEP 3: correlation multiple between substances and mixtures
dbcorrname = sprintf('dbcorrelation_mixture_%s.mat',o.mixturename);
if exist(fullfile(outputpath,dbcorrname),'file') 
    load(fullfile(outputpath,dbcorrname));
    if o.nmolchoose ~= size(dbcorrmix,1)
        dbcorrmix = nmrloadmatrix(moltestsortbyrho,mixturespec(:,2),nmrdb.dbpur,nmrdb.dbxpur,nmrdb.dbfit,'verbose',2,'maxlag',0.02);
        save(fullfile(outputpath,dbcorrname),'dbcorrmix')
    end
else
    dbcorrmix = nmrloadmatrix(moltestsortbyrho,mixturespec(:,2),nmrdb.dbpur,nmrdb.dbxpur,nmrdb.dbfit,'verbose',2,'maxlag',0.02);
    save(fullfile(outputpath,dbcorrname),'dbcorrmix')
end

%% STEP 4: inversion of matrix
prctiles = [0 25 40 50 60 75 85 90 95 100];
[moltest,imoltest] = intersect({nmrdb.dbcorrsub.dbcorrsubfull(:,1).subtestlist},moltestsortbyrho);
[~,isortbyrho] = intersect(moltestsortbyrho,moltest);
[~,imoltestsortbyrho] = intersect(isortbyrho,sort(isortbyrho));
imoltest = imoltest(imoltestsortbyrho);
dbcorrsubtest = nmrdb.dbcorrsub.dbcorrsubfull(imoltest,imoltest);
dbcorrmixtest = dbcorrmix(imoltestsortbyrho,1);

rang = nmrdb.dbclass.Tw(imoltest,nclasstheorical); % order of proximity of substances generated by cluster
splitsize = diff(find(diff([-inf;rang-(1:length(rang))'*0;inf]))); % find consecutive of number have a same groupe
classtmp = arrayfun(@(i) ones(1,splitsize((i)))*(i),1:length(splitsize),'UniformOutput',false);  class = [classtmp{:}]'; % conversion of rang to class which starts by class 1

proton = str2double(cellfun(@(m) uncell(regexp(regexp(nmrdb.dbpur.(m).formula,'H\d+','match'),'\d+','match')),moltestsortbyrho));
Mw = cellfun(@(m) sum(nmrdb.dbpur.(m).Mw),moltestsortbyrho);% molecualr mass of moltest
[~,~,dbout,dataplot] = nmrlsqnonneg(dbcorrsubtest,dbcorrmixtest,'nmoltest',o.nmolchoose,'proximity',class,'P',proton,...
                        'convert2molar','Mw',Mw,'dbfit',nmrdb.dbfit,'dbxpur',nmrdb.dbxpur,'dbpur',nmrdb.dbpur,'calibrationvalue',nmrdb.dbcalibration.calibcommon.b,...
                        'nmoltree',o.nmoltree,'prctiles',prctiles,'concmax',10000,'valuestoplot','conc_p50','fontsizegraph',9,'linewidthgraph',.5,'colormap',colormapgraph);  

%% output
R.dbout = dbout;
R.dataplot = dataplot;
R.rhotestsort = rhotestsort;
R.rang = rang;
R.mixturename = o.mixturename;
R.nplotfit = o.nplotfit;
R.corrmeanweightsort = corrmeanweightsort;
R.molsortbyrho = molsortbyrho;
R.dbinfo = dbinfo;

