function [dbmask,mask,dbmultiplet] = nmrloadmask(varargin)
%NMRLOADMASK load mask file in ods format 
% ods file containts sheets: function and its explication
%                            nmr multiplicity and its explication
%                            substance and its explication
%                            others worksheets correspond to the mask of substances reported in worksheet "substance"
% Syntax: db = nmrloadmask('parameter1',value1,'parameter2',value2,...)
%
% List of input pair parameters/values:
%         'path': full path where ODS file is located (default = find_path_toolbox('rmnspec'))
%      'odsfile': ODS file containing mask information (default='mask.ods')
% 'prefetchfile': file to be used for fast loading (default=dbmask.mat)
%
% OUTPUT
%           DBMASK: mask database structure containing
% dbmask.(mol).zone:	i{1...n}: included areas. e{1...n}: excluded areas
% dbmask.(mol).position:	chemical shift (ppm) of region (il corresponds to chemical shift of peak)
% dbmask.(mol).width1:	width of region
% dbmask.(mol).weight:	weight attributed to region (excluded zones = 0, included zones: 0-5)
% dbmask.(mol).buffer:	buffer in ppm to identify the beginning and the end of region
% dbmask.(mol).function:	chemical function corresponding to this region (peak)
% dbmask.(mol).proton:	number of proton corresponding to this region
% dbmask.(mol).multiplicity:	multiplet of peak
% dbmask.(mol).comment:	further information. comments
% dbmask.(mol).J1:	constant coupling 1
% dbmask.(mol).J2:	constant coupling 2
% dbmask.(mol).ppmmin:	Starting shift for the region
% dbmask.(mol).ppmmax:	Ending shift for the region
% dbmask.(mol).ppmwidth:	width of region (official)
% dbmask.(mol).gates: nx4 array, data for MULTIPB ([ppmmin ppmmax buffer weight])
%
%           MASK: all info contained in MASK.ODS file
%    DBMULTIPLET: ranking sub substances according to mulitplet types (s=singlet,d=doublet,t=triplet,q=quartet,p=pentet,st=sextet,dd=doublet of doublet
%                 dt=doublet of triplets,td=triplet of doublets,tt=triplet of triplets,mr=multiplet resolved,mnr=multiplet non resolved)
% dbmultiplet.(multipletrype) = {commonname1,commnname2}
%
% RMNSPEC v 0.1 - 23/01/13 - INRA\Olivier Vitrac, LNE\Mai Nguyen - rev. 27/03/13
% history
% 27/03/13 add output 'mask' and 'dbmultiplet'
% 01/11/13 ad field db.help
% default
default = struct('path',find_path_toolbox('rmnspec'),'odsfile','mask.ods','prefetchfile','dbmask.mat');
keyword = 'noprefetch';

% argcheck
o = argcheck(varargin,default,keyword);

% main
% check existance of database of mask 
if exist(fullfile(o.path,o.prefetchfile),'file') && ~o.noprefetch 
    dispf('MASK: reuse existing mask database')
    fileinfo(fullfile(o.path,o.prefetchfile))
    load(fullfile(o.path,o.prefetchfile)); % load database of mask
else % generate data
    mask = loadodsprefetch(fullfile(o.path,o.odsfile),'sheetname','all'); % load ods file of mask 
    molref = mask.substance.reference;
    molname = mask.substance.commonname;    
    maskfield = mask.nfo.fieldname(1:end-1); %exlude the last field = gates 
    nmaskfield = length(maskfield);
    maskprop = cell2struct(repmat({[]},nmaskfield,length(molref)),maskfield);
    dbmask = struct([]);
%     dbmaskhoriz = struct([]);
    for i = 1:length(molref)
        if isfield(mask,molref{i})==0; % if the mask of molecule is empty
           dispf('WARNING: the mask of ''%s'' is not available, please check',molref{i})
           dbmask(1).(molname{i}) = maskprop(i); 
        else
           for k = 1:nmaskfield
               maskprop(i).(maskfield{k}) = mask.(molref{i}).(maskfield{k});%dbmask(1).(molname{i}) = struct(maskfield{k},mask.(molref{i}).(maskfield{k}));               
           end  
           dbmask(1).(molname{i}) = maskprop(i);
           gates = zeros(length(dbmask.(molname{i}).weight),4);
           gates(:,1) = dbmask.(molname{i}).ppmmin;
           gates(:,2) = dbmask.(molname{i}).ppmmax;
           gates(:,3) = dbmask.(molname{i}).buffer;
           gates(:,4) = dbmask.(molname{i}).weight;
           dbmask.(molname{i}).gates = gates;
%            dbmaskhoriz(1).(molname{i}) = cell2struct(num2cell(struct2structtab(dbmask.(molname{i})),2),dbmask.(molname{i}).zone); 
        end 
    end
    % add help
    dbmask.help = cell2struct(mask.nfo.description,mask.nfo.fieldname);
    % save
    dispf('NMRLOADMASK:\t new/updated prefetch file')
    save(fullfile(o.path,o.prefetchfile),'dbmask','mask') %,'dbmaskhoriz'
    fileinfo(fullfile(o.path,o.prefetchfile))
end
% Additional output: dbmultiplet (indexing according to each multiplet type
if nargout>2
    fdb = fieldnames(rmfield(dbmask,'help'));
    multiplet = cellfun(@(m) dbmask.(m).multiplicity,fdb,'UniformOutput',false);
    multiplettype = mask.multiplicity.multiplicity;
    dbmultiplet = cell2struct(repmat({[]},length(multiplettype),1),multiplettype);
    for i = 1:length(multiplettype)
        substancetostore = cell(length(fdb),1);
        for j = 1:length(fdb)
            clear any, clear tmp
            tmp = cellfun(@(x) ischar(x) && strcmp(x,multiplettype{i,1}), multiplet{j,1});
            if any(tmp > 0); substancetostore{j,1} = fdb(j); end
        end
        substancetostore(cellfun(@isempty,substancetostore)) = [];
        dbmultiplet.(multiplettype{i}) = substancetostore;
    end
end