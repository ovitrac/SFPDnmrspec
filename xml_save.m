function xml_save( file, S, att_switch, header )
%function xml_save( filename, S, att_switch )
%
% Saves structure or variable(s) S in xml format to a file.
%
% INPUT 
%   filename     filename 
%   S            Matlab variable or structure to store in file.
%   att_switch   optional, 'on' stores XML type attributes (default), 
%                'off' doesn't store XML type attributes
%   header       string or cell array of strings to set the header
%               default value = sprintf('SFPP3\\INRA-Olivier Vitrac %s', datestr(now))
%
% OUTPUT 
%   none
%
% RELATED 
%   xml_load, xml_format, xml_parse, (xmlread, xmlwrite)
%
% --------------------------------------------------------------
% Copyright Geodise/GEM 2002-2003, University of Southampton, UK
% Author: Marc Molinari <m.molinari@soton.ac.uk>
% $Revision: 2.0.1.1 $ $Date: 2003/12/11 14:26:32 $

% Revised INRA\Olivier Vitrac to handle n-D cell arrays and structure arrays, better extension control, multiple line headers
% $$ Last revision by INRA\Olivier Vitrac 27/11/2014 $$

% Please read the licence in file licence.txt.

%---------------------------
% DEFAULT
att_switch_default = 'on';
header_default = sprintf('SFPP3\\INRA-Olivier Vitrac %s', datestr(now));

% INIT 
if (nargin<2) || ~ischar(file)
  error('%s requires 2 or 3 parameters: filename and variable, optionally att_switch.',mfilename);
end

if nargin<3 ,att_switch=att_switch_default; end
if nargin<4, header = ''; end
if isempty(header), header = header_default; end
if ~iscell(header), header = {header}; end
nheader = length(header);
if isempty(att_switch), att_switch=att_switch_default; end
if ~strcmpi(att_switch, 'off') && ~strcmpi(att_switch, 'idx'); att_switch = 'on'; end

%-----------------------------------------------
% append '.xml'
[~,~,xmlext] = fileparts(file);
if ~strcmpi(xmlext,'.xml'), file = strcat(file, '.xml'); end

%-----------------------------------------------
% overwrite ?
%if (exist(file))
%  q = sprintf('Overwrite existing file %s ?', file);
%  title = 'Warning: File already exists!';
%  dosave = questdlg(q, title, 'No', 'Yes', 'No');
%  if strcmpi(dosave, 'no')
%    return;
%  end
%end

%=====================================================
fid = fopen(file, 'w');
if fid==-1
  error(['Error while writing file ', file]);
end

% write file header
fprintf(fid, '<?xml version="1.0"?>\n');
% fprintf(fid, sprintf('<!-- SFPP3\\INRA-Olivier Vitrac %s -->\n', datestr(now)));
% fprintf(fid, '<!-- SFPP3\\INRA-Olivier Vitrac %s -->\n', datestr(now));
for i = 1:nheader 
    if i==1, fprintf(fid,'<!-- '); end
    fprintf(fid,'%s',header{i});
    if (nheader>1 && i<nheader),  fprintf(fid,'\n'); end
    if i==nheader, fprintf(fid,' -->\n'); end
end

% write XML string
fprintf(fid, '%s', xml_format(S, att_switch));

% close file
fclose(fid);
