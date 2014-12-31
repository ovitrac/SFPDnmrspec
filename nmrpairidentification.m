function pairidentification = nmrpairidentification(nmrdb,mixturespec,varargin)
%NMRPAIRIDENTIFICATION calculates pairwise correlation coefficient between 1H NMR mixture spectrum and a dictionary of 1H NMR spectra
%   R = nmrpairidentification(nmrdb,mixture,'property1',value1,'property2',value2,...)
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
%
% OUTPUT
%   R = structure with fields 
%     .corrmeanweightsort: correlation coefficient averaged, weighted and sorted in descending order
%     .molsortbyrho: substance names according to rank of corr coeff

%
% See also: deformulateNMR, nmrlsqnonneg, nmrloadbase, nmrloadascii, nmrloaddbfit, nmrloadmatrix, nmrloadcorrcoeff, plotnmrlsqnonneg

% RMNSPEC v 0.5 - 06/11/2014 - INRA\Olivier Vitrac, LNE\Mai Nguyen - rev. 26/11/14
% History
% 26/11/14: change output structure and change extension of dbname

% default
outputpathdefault = pwd;
default = struct('worse',10,...
                 'idsubstance','commonname',... 
                 'threshold',0.6,... 
                 'outputpath',outputpathdefault,...
                 'mixturename','mixtureS1',...
                 'dbnameext','');

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
% extract substance names for calculation : commonname
mol = fieldnames(rmfield(nmrdb.dbfit,'help'));
repet = cellfun(@isempty,regexp(mol,'R\d+'));% extract only fitmol without repetition
mol = mol(repet); nmol = length(mol);
% output directory
outputpath = fullfile(o.outputpath,o.mixturename);
if ~exist(outputpath,'dir'), mkdir(outputpath); end 
% calculation of pairwise correlation coefficient 
if isempty(o.dbnameext)
    dbcorrcoeffname = sprintf('%s.identification.mat',o.mixturename);
else
    dbcorrcoeffname = sprintf('%s%s',o.mixturename,o.dbnameext);
end
if exist(fullfile(outputpath,dbcorrcoeffname),'file')
    dispf('\tload data of pair identification'), fileinfo(fullfile(outputpath,dbcorrcoeffname))
    pairidentification = load(fullfile(outputpath,dbcorrcoeffname));
else
    db = nmrloadcorrcoeff(mol,mixturespec(:,2),nmrdb.dbpur,nmrdb.dbxpur,nmrdb.dbfit,'worse',o.worse,'idsubstance',o.idsubstance,'threshold',o.threshold);
    % ranking
    corrmeanw = cell2mat(arrayfun(@(m) db.corrmeanweight{(m),1}(1,3),1:nmol,'UniformOutput',false))';
    [corrmeanweightsort,iwsort] = sort(corrmeanw,'descend'); %#ok<UDIM>
    molsortbyrho = db.subtestlist(iwsort);
    % save
    if ~exist(outputpath,'dir'), mkdir(outputpath); end
    pairidentification.corrmeanweightsort = corrmeanweightsort;
    pairidentification.molsortbyrho = molsortbyrho;
    pairidentification.mixturename = o.mixturename;
    pairidentification.dbinfo = dbinfo;
    save(fullfile(outputpath,dbcorrcoeffname),'-struct','pairidentification')
end     
