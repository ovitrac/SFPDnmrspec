function xml_help
% --------------------------
% XML TOOLBOX FOR MATLAB 2.0
% --------------------------
%
% FUNCTIONS:
%  xml_format       converts a Matlab variable/structure into XML string
%  xml_parse        parses & converts XML string into Matlab variable/structure
%  xml_save         saves a Matlab variable/structure in XML format in a file
%  xml_load         loads the .xml file written with xml_save back into a variable
%  xml_help         this file, displays info about available xml_* commands
%  tests/xml_tests  tests the xml toolbox by writing/reading a number of xml test files
%  strsplit         utility function which splits a string at specified characters
%
% FILES:
%  doc/xml_matlab   documentation containing info on installation, usage, implementation, etc.
%  matlab.xsd       contains the Schema to validate xml files for the toolbox 
%                   (if not present, look at http://www.geodise.org/matlab.xsd)
%
% RELATED:
%  xmlread, xmlwrite (shipped with Matlab from version 6.5)
%
% Further information can be obtained by using the help command on
% a specific function, e.g. help xml_format.
%
% --------------------------------------------------------------
% Copyright Geodise/GEM 2002-2003, University of Southampton, UK
% Author: Marc Molinari <m.molinari@soton.ac.uk>
% $Revision: 2.0.1.1 $ $Date: 2003/12/11 14:26:32 $

% Please read the licence in file licence.txt.

% VERSION CONTROL
%  0.1, 04/12/2002 created by <m.molinari@soton.ac.uk> (mm)
%  1.0, 20/12/2002 Release 1.0
%  2.0, 20/10/2003 Release 2.0

help xml_help
