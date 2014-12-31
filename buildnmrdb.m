function buildnmrdb(loadingfile,massfile,varargin)
%BUILDNMRDB builds nmrdb for opensource project from mass ans loading files
% buildnmrdb(massfile,loadingfile,'property1',value1,'property2',value2,...)
% INPUTS
%      massfile: fullpath and file of mass matrix (supplied in .mat)
%   loadingfile: fullpath and file of loading base (supplied in .xml) 
%
% Optional pair property/value (default value)
%  'outputpath': location to save nmrdb (default = pwd)
%  'dbfilename': name file of db (default = 'nmrdb.mat')
%
% RMNSPEC v 0.5 - 07/11/2014 - INRA\Mai Nguyen, Olivier Vitrac  - rev. 27/11/14
% History
% 27/11/14: add control on mass and loading files

% default
outputpathdefault = pwd;
default = struct('outputpath',outputpathdefault,'dbfilename','nmrdb.mat');
             
% ARGCHECK
o = argcheck(varargin,default);
if ~exist(loadingfile,'file'), error('The supplied loading file "%s" does not exist',loadingfile), end
if ~exist(massfile,'file'), error('The supplied loading file "%s" does not exist',massfile), end
% reconstitute NMRDB
if ~exist(fullfile(o.outputpath,o.dbfilename),'file')
    % Load loading base suplied in xml format: dbpur, dbxpur, dbfit, dbclass
    nmrdb = xml_load(loadingfile);
    % control
    if ~isstruct(nmrdb) || ~isfield(nmrdb,'dbpur') || ~isfield(nmrdb,'dbxpur') || ~isfield(nmrdb,'dbfit') || ~isfield(nmrdb,'dbclass') || ~isfield(nmrdb,'dbcalibration') || ~isfield(nmrdb,'dbcorrsub')
        error('The supplied loading file is invalid. Please contact administrator')
    end
    mol = fieldnames(rmfield(nmrdb.dbpur,'help')); nmol = length(mol);
    molsort = sort(mol);
    % reconstitute ppm vector
    nmrdb.dbxpur.ppm = linspace(nmrdb.dbxpur.ppm(1),nmrdb.dbxpur.ppm(2),nmrdb.dbxpur.ppm(3))';
    % reconstitute signals in dbfit
    gaussiankernel = @(x,position,width)  exp(-((x-position)./(0.6006.*width)) .^2);
    k = 1;
    for imol = 1:nmol
        for j = 1:size(nmrdb.dbfit.(mol{imol}).all,1)
            x = nmrdb.dbfit.(mol{imol}).all(j).ppm(1):nmrdb.dbfit.(mol{imol}).all(j).ppm(3):nmrdb.dbfit.(mol{imol}).all(j).ppm(2);
            position = nmrdb.dbfit.(mol{imol}).all(j).position;
            width =  nmrdb.dbfit.(mol{imol}).all(j).width;
            weight =  nmrdb.dbfit.(mol{imol}).all(j).weight';
            n = nmrdb.dbfit.(mol{imol}).all(j).npeaktofit;
            nmrdb.dbfit.(mol{imol}).all(j).Ifit = gaussiankernel(repmat(x(:),1,n),position(k*ones(1,numel(x)),1:n),width(k*ones(1,numel(x)),1:n))*weight(1:n,k);
            nmrdb.dbfit.(mol{imol}).all(j).ppm = x;
        end
    end
    
    % load mass matrix supplied in .mat
    load(massfile)
    % reconstitute mass matrix -> structure dbcorrsub.dbcorrsubfull with fields 'corrmax' 'corrzerolag' 'rho2' 'subreflist' 'subtestlist'    
    corrmax = cell(nmol,nmol);
    k = 1;
    for i = 1:nmol
        npeak = length(dbcorrsubtable(dbcorrsubtable(:,1)==i,1)); 
        for j = 1:nmol, corrmax{i,j} = dbcorrsubtable(k:k+npeak-1,j+1); end
        k = k+npeak;
    end
    dbcorrsubfull = repmat(struct('corrmax',[],'corrzerolag',[],'rho2',[],'subreflist','','subtestlist',''),nmol,nmol);
    for i = 1:nmol
        for j = 1:nmol
            dbcorrsubfull(i,j).corrmax = corrmax{i,j};
            dbcorrsubfull(i,j).subreflist = molsort{j};
            dbcorrsubfull(i,j).subtestlist = molsort{i};
        end
    end
    nmrdb.dbcorrsub.dbcorrsubfull = dbcorrsubfull;
    % save nmrdb in .mat format in tempdir 
    save(fullfile(o.outputpath,o.dbfilename),'nmrdb')
end
