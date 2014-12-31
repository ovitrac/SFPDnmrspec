function [X] = xml_load( file )
%function [X] = xml_load( file )
%
% Loads XML file and converts it into Matlab structure or variable.
%
% INPUT 
%   file     filename of xml file written with xml_save
%
% OUTPUT 
%   X        structure variable containing file contents
%  [xmlV     !changed: returns empty string to maintain compatibility]
%
% SEE ALSO 
%   xml_format, xml_parse, xml_save, (xmlread, xmlwrite)
%
% --------------------------------------------------------------
% Copyright Geodise/GEM 2002-2003, University of Southampton, UK
% Author: Marc Molinari <m.molinari@soton.ac.uk>
% $Revision: 2.0.1.1 $ $Date: 2003/12/11 14:26:32 $

% $$ Last revision by INRA\Olivier Vitrac 10/10/2014 $$

% Please read the licence in file licence.txt.

% INTERFACE CHANGE:
% xml_load does not return XML version as second output any more.

% ----------------------------------------------------------
% Initialisation and checks on input parameters

% set XML TB version number
xml_tb_version = '2.01';

% check input parameters
if (nargin<1)
  error([mfilename,' requires 1 parameter: filename.']);
end

% append '.xml'
[~,~,xmlext] = fileparts(file);
if ~strcmpi(xmlext,'.xml'), file = strcat(file, '.xml'); end

%-----------------------------------------------
% check existence of file
if ~exist(file,'file')
  error([mfilename, ': could not find ', file]);
end

%-----------------------------------------------
fid = fopen(file, 'r');
if fid==-1
  error(['Error while opening file ', file, ' for reading.']);
end

% parse file content into blocks
str = char( fread(fid)' ); % read in whole file
fclose( fid );

if (length(str)<3)
  error([file, ' does not seem to be a valid xml file.']);
end

%-----------------------------------------------
% parse content, identify blocks
X = xml_parse(str);

return
