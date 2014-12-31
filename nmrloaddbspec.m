function varargout = nmrloaddbspec(dblist,varargin)
%NMRLOADDBSPEC loads chosen nmr spectral databases
% SYNTAX [db1 db2] = nmrloaddbspec({'db1' 'db2'})
% INPUTS
%   dblist: string, list of name of database to be loaded (if empty = 'all')
% VARARGIN
% 'dbpath': direction to store full database (.mat)
% 'dbname': string, name of full database (.mat) containing spectral databases (default='nmrdatabase.mat')
% OUTPUT
% ----------- HELP TO DO--------------
%
% RMNSPEC v 0.1 - 04/11/13 - INRA\Olivier Vitrac, LNE\Mai Nguyen - rev
% history

% default
default = struct('dbpath',find_path_toolbox('rmnspec'),'dbname','nmrdatabase.mat');
% argcheck
o = argcheck(varargin,default);
if nargin<1, dblist = {}; end

% check name of db
if ~isempty(dblist) && ~iscell(dblist), dblist = {dblist}; end
dbinside = who('-file',fullfile(o.dbpath,o.dbname));
if ~isempty(dblist), ndb = length(dblist); 
else ndb = length(dbinside);
end
[~,ifound] = intersect(dbinside,dblist);
if length(ifound)<length(dblist)
    dispf('ERROR\t%d database(s) cannot be found in %s',length(ifound)-ndb,fullfile(o.dbpath,o.dbname))
    cellfun(@(m) dispf('\t''%s'' is missing',m),setdiff(dblist,{dbinside}))
    error('%d missing values (see above)',length(ifound)-ndb)
end                            

% load db
varargout = cell(ndb,1);
if isempty(dblist)
    for i = 1:ndb, varargout{i} = load(fullfile(o.dbpath,o.dbname),dbinside{i}); end
else
    for i = 1:ndb, varargout{i} = load(fullfile(o.dbpath,o.dbname),dblist{i}); end
end
    

    