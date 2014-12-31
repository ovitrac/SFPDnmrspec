function tmp = loadSFPD(databasefile)
%LOADSFPP3 load a database in XML (force XML use even if .mat requested)
%   example: tmpSFPP3 = loadSFPP3('QSPR_arc2006.mat'); [arc,alive]=deal(tmpSFPP3.arc,tmpSFPP3.alive);

% MS-MATLAB-WEB 3.0 - 26/01/09 - Olivier Vitrac - rev.

% definitions
xml = '.xml';

%arg check
if nargin~=1, error('one argument is required'), end

% Conversion (if required) and load
[p,n,e] = fileparts(databasefile);

if isempty(e) || strcmpi(e,'.mat') % MAT asked but .XML prefered in SFPP3
    tmpfile = fullfile(p,[n xml]);
    if ~exist(tmpfile,'file') % try a conversion in XML
        tmp = load(databasefile);
        xml_save(tmpfile,tmp);
    else
        tmp = xml_load(tmpfile);
    end
elseif strcmpi(e,xml) && exist(databasefile,'file')
    tmp = xml_load(databasefile);
else
    error('LOADSFPD: invalid file ''%s''',databasefile)
end
