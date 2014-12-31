function R = nmrdeconvolution(nmrdb,mixturespec,pairidentification,varargin)
%NMRDECONVOLUTION deconvolves a 1H NMR mixture spectrum according to a dictionary of 1H NMR spectra from ranking of pairwise correlation coefficient
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
%       'outputpath': string, directory of output for saving results database 
%      'mixturename': string, name of mixture, for legend (default = 'mixtureS1')
%       'nmolchoose': numeric, number of substances choosen for deconvolution, higher precedence than 'rhothreshold' (default=[])
%     'rhothreshold': numeric < 1, threshold of pair correlation coefficient for choosing number of substances for deconvolution (default=0.75)
%  'nclasstheorical': numeric, number of classes in spectral database of substances (default=[])
%'rhothresholdclass': numeric, threshold of corr.coeff to sort substances by class (default=0.75)
%         'nmoltree': number substances to be plotted in scenario tree (default=5)
%         'nplotfit': numeric, number of substances to be fitted and ploted
%    'colormapgraph': colormap for tree graph (default=jet(64) in log scale)
%
% OUTPUT
%   R = structure with fields
%
% PRINCIPLE OF THE METHOD
% List of tasks carried out in DEFORMUALTENMR
% step 1: extract the best "10" substances or rho > 0.75 or ... 
% step 2: multiple correlation (nmrloadmatrix)
% step 3: inversion of matrix (nmrlsqnonngeg)
%
% See also: nmrpairidentification, deformulateNMR, nmrlsqnonneg, nmrloadbase, nmrloadascii, nmrloaddbfit, nmrloadmatrix, nmrloadcorrcoeff, plotnmrlsqnonneg

% RMNSPEC v 0.5 -06/11/2014 - INRA\Olivier Vitrac, LNE\Mai Nguyen - rev. 26/11/14
% History
% 26/11/14: change output structure and change dbname

% default
outputpathdefault = pwd;
default = struct('outputpath',outputpathdefault,...
                 'mixturename','mixtureS1',...
                 'nmolchoose',[],...
                 'rhothreshold',0.75,...
                 'dbnameext','',...
                 'nclasstheorical',[],...
                 'rhothresholdclass',0.75,...
                 'nmoltree',5,...
                 'nplotfit',4,...
                 'concmax',10000,...
                 'colormapgraph',[]);
% ARGCHECK
o = argcheck(varargin,default,'nostructexpand','case');
% MANDATORY PROPERTIES: nmrdb,mixture
if nargin < 3, error('3 arguments are required'), end
if ischar(nmrdb)
    dbinfo = fileinfo(nmrdb);
    if exist(nmrdb,'file')
        dispf('\tload reference NMR database:'), fileinfo(nmrdb), load(nmrdb)
    else
        error('the reference NMR database ''%s'' does not exist in ''%s''',lastdir(nmrdb),rootdir(nmrdb))
    end
end
if ~isstruct(nmrdb) || ~isfield(nmrdb,'dbpur') || ~isfield(nmrdb,'dbxpur') || ~isfield(nmrdb,'dbfit') || ~isfield(nmrdb,'dbclass') || ~isfield(nmrdb,'dbcalibration') || ~isfield(nmrdb,'dbcorrsub')
    error('nmrdb must be created with buildnmrdb function')
end
if ischar(mixturespec)
    if exist(mixturespec,'file')
        dispf('\tload mixture spectrum:'), fileinfo(mixturespec), mixturespec = load(mixturespec);
    else
        error('the mixture spectum ''%s'' does not exist in ''%s''',lastdir(mixturespec),rootdir(mixturespec))
    end
end
if ~isnumeric(mixturespec) || size(mixturespec,2)~=2, error('mixture must be nx2 numeric vector (1st column: ppm and 2nd column: I'), end
if ~isnumeric(mixturespec) || size(mixturespec,1)~=length(nmrdb.dbxpur.ppm) || mixturespec(1,1)~=nmrdb.dbxpur.ppm(1) || mixturespec(end,1)~=nmrdb.dbxpur.ppm(end)
    error('Please define mixture NMR spectrum with chemical shift limited in %d-%d ppm and containing %d points',nmrdb.dbxpur.ppm(1),nmrdb.dbxpur.ppm(end),length(nmrdb.dbxpur.ppm))
end
if ischar(pairidentification)
    if exist(pairidentification,'file')
        dispf('\tload data of pair identification'), fileinfo(pairidentification)
        pairidentification = load(pairidentification);
    else
        error('datafile of pair identification ''%s'' does not exist in ''%s''',lastdir(pairidentification),rootdir(pairidentification))
    end
end
if ~isstruct(pairidentification) || ~isfield(pairidentification,'corrmeanweightsort') || ~isfield(pairidentification,'molsortbyrho') || ~isfield(pairidentification,'mixturename') 
    error('''pairidentification'' must be created with NMRPAIRIDENTIFICATION function')
end
 
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
    nclasstheorical = ceil(csaps(nmrdb.dbclass.rhowfilt(:,3),1:size(nmrdb.dbclass.rhowfilt,1),[],o.rhothresholdclass));
else warning('No information about the classification of spectral database. Number of classes il the number of substances: %d',nmol)
end
% output directory
outputpath = fullfile(o.outputpath,o.mixturename);
if ~exist(outputpath,'dir'), mkdir(outputpath); end 

%% calculation
if isempty(o.dbnameext)
    dbcorrname = sprintf('%s.deconvolution.mat',o.mixturename);
else
    dbcorrname = sprintf('%s%s',o.mixturename,o.dbnameext);
end

if exist(fullfile(outputpath,dbcorrname),'file') 
    dispf('\tload data of deconvolution:'), fileinfo(fullfile(outputpath,dbcorrname)), 
    R = load(fullfile(outputpath,dbcorrname));
else
    % STEP 1: choose of substances for deconvolution 
    molsort = pairidentification.molsortbyrho ;
    rhosort = pairidentification.corrmeanweightsort;
    if ~isempty(o.rhothreshold)
        tmpmol = molsort(1:length(rhosort(rhosort>o.rhothreshold)));
        tmpmol = setdiff(tmpmol,{'PE' 'PP' 'PS'}');
    end
    if ~isempty(o.nmolchoose)
        tmpmol = molsort(1:o.nmolchoose);
        tmpmol = setdiff(tmpmol,{'PE' 'PP' 'PS'}');
        lengthtest = o.nmolchoose+o.nmolchoose-length(tmpmol);
        while length(tmpmol) < o.nmolchoose
            tmpmol = molsort(1:lengthtest);
            tmpmol = setdiff(tmpmol,{'PE' 'PP' 'PS'}');
            lengthtest = lengthtest+1;
        end 
    end
    moltestsortbyrho = molsort(ismember(molsort,tmpmol)); % list of n mols test in decreasing order of rho2
    rhotestsort = rhosort((ismember(molsort,tmpmol)),1);% list of rho2 of n mols test

    % STEP 2: correlation multiple between substances and mixtures
    dbcorrmix = nmrloadmatrix(moltestsortbyrho,mixturespec(:,2),nmrdb.dbpur,nmrdb.dbxpur,nmrdb.dbfit,'verbose',2,'maxlag',0.02);
        
    % STEP 3: inversion of matrix
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
                            'nmoltree',o.nmoltree,'prctiles',prctiles,'concmax',o.concmax,'valuestoplot','conc_p50','fontsizegraph',9,'linewidthgraph',.5,'colormap',colormapgraph);  

    % output
    R.dbcorrmix = dbcorrmix;
    R.dbout = dbout;
    R.dataplot = dataplot;
    R.rhotestsort = rhotestsort;
    R.rang = rang;
    R.mixturename = o.mixturename;
    R.nplotfit = o.nplotfit;
    R.dbinfo = dbinfo;
    % save
    if ~exist(outputpath,'dir'), mkdir(outputpath); end
    save(fullfile(outputpath,dbcorrname),'-struct','R')
end

