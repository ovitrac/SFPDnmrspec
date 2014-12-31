function saveSFPD(databasefile,vars,varnames)
%SAVESFPP3 save a database in XML (even if a '.mat' is requested)
%   example: saveSFPD('nmrdb.mat',{dbpur dbfit},{'dbpur' 'dbfit'})

% MS-MATLAB-WEB 3.0 - 26/01/09 - Olivier Vitrac - rev.

% definitions
xml = '.xml';

%arg check
if nargin~=3, error('three arguments are required'), end
if ~iscell(vars), vars = {vars}; end
[p,n,e] = fileparts(databasefile);
if ~strcmpi(e,xml);
    tmpfile = fullfile(p,[n xml]);
else
    tmpfile = databasefile;
end    

%save
dim = find(size(vars)>1,1,'first');
if isempty(dim), dim = 2; end
tmp = cell2struct(vars,varnames,dim);
xml_save(tmpfile,tmp);