function [X, xmlV] = xml_parse(str, att_switch, X, level)
%function [X, xmlV] = xml_parse(str, [att_switch], [X], [level])
%
% Parses XML string str and returns matlab variable/structure.
% This is a non-validating parser!
%
% INPUT 
%   str         xml string, possibly from file with function xml_load.m
%   att_switch  'on'- reads attributes, 'off'- ignores attributes
%   X           Optional. Variable which gets extended or whose
%               substructure parameters get overridden by entries in
%               the string.
%   level       Internal xml level. Should not be used by user.
%
% OUTPUT 
%   X        Matlab variable or structure 
%  [xmlV     not supported any more, contained the xml version definition]
%
% RELATED 
%   xml_format, xml_load, xml_save, (xmlread, xmlwrite)
%
% --------------------------------------------------------------
% Copyright Geodise/GEM 2002-2003, University of Southampton, UK
% Author: Marc Molinari <m.molinari@soton.ac.uk>
% $Revision: 2.0.1.1 $ $Date: 2003/12/11 14:26:32 $

% INRA\Olivier Vitrac, add UTF-8 decoding, add single and all integer types
% $$ Last revision by INRA\Olivier Vitrac 10/10/2014 $$

% Please read the licence in file licence.txt.

% ----------------------------------------------------------
% Initialisation and checks on input parameters

% set XML TB version number
xml_tb_version = '2.01';

% we need function name for subsequent recursive calls
thisfunc = mfilename;

% define persistent (static) variable in which read-in
% xml_tb_version number gets stored
persistent xmlTBVersion;

% check input parameters
if (nargin<1), error([mfilename, ' needs at least 1 input parameter.']); end
% default attribute switch setting
if ((nargin<2) || ~strcmp(att_switch, 'off')), att_switch = 'on'; end
% default base variable setting
if (nargin<3),  X = struct([]); end
% default root level setting
if (nargin<4), level = 0; end
% check input string
if isempty(str), return; end

% define variables
xmlVersion = '';

%---------------------------
% remove all <! execute and comment entries from str by blanking out
execpos = strfind2(str, '<!');
if ~isempty(execpos)
  allclose = strfind2(str, '>');
  for x=1:length(execpos)
    xstart   = execpos(x);
    idxclose = find(allclose > xstart);
    xend     = allclose(idxclose(1));
    str(xstart:xend) = blanks(xend-xstart+1);
  end
end

%---------------------------
% find xml string elements
tmp = strfind2(str, '<');

pclose = sort( [strfind2(str, '</'), ...
    strfind2(str, '/>'), ...
    strfind2(str, '-->'), ...
    strfind2(str, '?>')] );
popen  = setdiff(tmp, pclose);

% check for correct number of start and end tags
if length(popen) ~= length(pclose)
  error('XML parse error: Number of element start and end tags does not match.');
end

np = length(popen);
pidx = sortrows([popen(:), -ones(np,1); ...
    pclose(:), ones(np,1)]);

% loop through all elements identified (on level 0 only root 
% which will call further instances of this function)
i=1;
sumparenths = 0;

while i<=size(pidx,1)
  
  entrystart = pidx(i,1);
  sumparenths = sumparenths + pidx(i,2);
  while sumparenths ~= 0
    i = i+1;
    sumparenths = sumparenths + pidx(i,2);
  end
  entryend = pidx(i,1);
  tmp = str(entrystart+1:entryend-1);
  
  TYPE = ''; NAME = ''; IDX=[]; SIZE=[0 0]; FIELDS=[]; TAG = ''; %#ok<*NASGU>
  
  headsep = strfind2(tmp, '>');
  if isempty(headsep)
    % deal with "/>" empty elements by using the whole tmp string
    headsep = length(tmp);
  end
  
  namesep = min([strfind2(tmp, ' '), strfind2(tmp, '>')]);
  if isempty(namesep)
    TAG = tmp;
  else
    TAG = tmp(1:namesep-1);  
  end
  
  header  = tmp(namesep+1:headsep-1);
  content = tmp(headsep+1:end);
  
  % make sure that we have size [0 0] and not [1 0]
  if isempty(content)
    content = '';
  end
  
  % parse header for attributes
  att_lst = header;
  att_lst(att_lst=='=') = '"';
  h = strsplit(att_lst, '"');
  
  if strcmp(att_switch, 'on')
    for k=1:2:length(h)-1
      switch(h{k})
        case 'idx'
          IDX = str2double(h{k+1});
        case 'name'
          NAME = h{k+1};
        case 'size'
          SIZE = str2num(h{k+1});
        case 'fields'
          FIELDS = strsplit(h{k+1},' ');
        case 'type'
          TYPE = h{k+1};
        %case 'value'
          % deal with this case like with any other attributes
          % in XML Toolbox Version 3.0
        case 'xml_tb_version'
          xmlTBVersion = str2double(h{k+1});
      end
    end
  end
  
  % special names
  switch (TAG(1))
    case {'?', '!'}   
      % ignore entity declarations and processing instructions
      % Note: we also ignore the <?xml ...> entry with version number.
      i=i+1;
      continue;
  end
  
  if isempty(xmlTBVersion) && (level==0)
    % this is possibly a version 1.x XML string
    if (strcmp(TAG, 'struct')  || ...
        strcmp(TAG, 'double')  || ...
        strcmp(TAG, 'char')    || ...
        strcmp(TAG, 'boolean') || ...
        strcmp(TAG, 'complex') || ...
        strcmp(TAG, 'sparse')  || ...
        strcmp(TAG, 'cell')    || ...
        strcmp(TAG, 'single')  || ...
        strcmp(TAG, 'uint8')   || ...
        strcmp(TAG, 'uint16')  || ...
        strcmp(TAG, 'uint32')  || ...
        strcmp(TAG, 'uint64')  || ...
        strcmp(TAG, 'int8')    || ...
        strcmp(TAG, 'int16')   || ...
        strcmp(TAG, 'int32')   || ...
        strcmp(TAG, 'int64')     ... 
       )
      xmlTBVersion = 1.0;
      NAME = 'root';
    else
      % att_switch is probably set to 'off'
      xmlTBVersion = 2.0;
    end
  end
  
  if (xmlTBVersion >= 2.0)
    % from version 2.0 we have NAME = TAG and TYPE is 
    % usually given, except when using att_switch=off
    NAME = TAG;
    if isempty(TYPE)
      TYPE = 'char';
    end
  else % (xmlTBVersion < 2.0)
    % version 1.0 has type as tag and name as attribute.
    % if no name is given, assign 'item'
    TYPE = TAG;
    if isempty(NAME)
      NAME = 'item';
    end
  end
  
  % remove namespace from NAME
  f = strfind2(NAME, ':');
  if ~isempty(f)
    NAME = NAME(f+1:end);
  end
  
  % remove namespace from TYPE
  f = strfind2(TYPE, ':');
  if ~isempty(f)
    TYPE = TYPE(f+1:end);
  end
  
  % make sure TYPE is valid
  if isempty(NAME) || isempty(TYPE)
    error('NAME or TYPE is empty!')
  end
  
  % check if type is correct
  if strcmp(TYPE, 'char') && any(strfind2(content, '<'))
    if strcmp(att_switch, 'on')
      TYPE = 'struct';
    else
      TYPE = 'parent';
    end
  end
  
  % check if index is correct
  if IDX==0
    IDX = [];
  end
  
  if ~isempty(X) && isfield(X, NAME) && isempty(IDX)
    cont_list = {X.(NAME)};
    found = 0;
    % this loop makes sure that the current entry is inserted 
    % after the last non-empty entry in the content vector cont_list
    for cc=length(cont_list):-1:1
      if ~isempty(cont_list{cc})
        found=1;
        break
      end
    end
    if ~found
      IDX = max(cc-1,1);
    else
      IDX = cc+1;
    end
  end
  
  if isempty(IDX) && ~isempty(X) && strcmp(NAME, 'item')
    % make sure that when we have a character array the IDX of the
    % new vector is set to 2 and not to the end+1 index of the string.
    if isa(X, 'char')
      IDX = 2;
    else
      IDX = length(X)+1;
    end
  end
  
  if isempty(IDX)
    % if everything else did not produce a result, assign IDX=1
    IDX = 1;
  end
  
  % switch board which decides how to convert contents according to TYPE
  switch lower(TYPE)
    
    % ========================
    case '?xml'
      % xml version definition
      % NOTE: this is never reached from xml_tb_version >= 2.0
      xmlVersion = content;
      
      % ========================
    case '!--'
      % comment, just ignore
      i = i+1;
      continue
      
      % ========================
    case {'double' 'integer', 'float', 'numeric'}
      c = str2num(content); if ~prod(SIZE)==0, c = reshape(c, SIZE); end %#ok<*ST2NM>
    case 'single'
      c = single(str2num(content)); if ~prod(SIZE)==0, c = reshape(c, SIZE); end
    case 'uint32'
        c = uint32(str2num(content)); if ~prod(SIZE)==0, c = reshape(c, SIZE); end
    case 'uint16'
        c = uint16(str2num(content)); if ~prod(SIZE)==0, c = reshape(c, SIZE); end
    case 'uint8'
        c = uint8(str2num(content)); if ~prod(SIZE)==0, c = reshape(c, SIZE); end
    case 'uint64'
        c = uint64(str2num(content)); if ~prod(SIZE)==0, c = reshape(c, SIZE); end
    case 'int32'
        c = int32(str2num(content)); if ~prod(SIZE)==0, c = reshape(c, SIZE); end
    case 'int16'
        c = int16(str2num(content)); if ~prod(SIZE)==0, c = reshape(c, SIZE); end
    case 'int8'
        c = int8(str2num(content)); if ~prod(SIZE)==0, c = reshape(c, SIZE); end
    case 'int64'
        c = int64(str2num(content)); if ~prod(SIZE)==0, c = reshape(c, SIZE); end
      
      % ========================
    case {'logical', 'boolean'}
      c = logical(str2num(content));
      if ~prod(SIZE)==0, c = reshape(c, SIZE); end
      
      % ========================  
    case {'char', 'string'}
      c = charunsubst(content);
      if isempty(c) && (length(c) ~= prod(SIZE))
        % this is a string containing only spaces
        c = blanks(prod(SIZE));
      end
      
      if ~prod(SIZE)==0
          [utf8start,utf8stop] = regexp(c,'&#\d*;');
          if any(utf8start)
              for iutf8=1:length(utf8start)
                  c = [c(1:utf8start(iutf8)-1) char(str2double(c(utf8start(iutf8)+2:utf8stop(iutf8)-1))) c(utf8stop(iutf8)+1:end)];
                  nutf8 = (utf8stop(iutf8)-utf8start(iutf8));
                  utf8start = utf8start  - nutf8;
                  utf8stop  = utf8stop   - nutf8;
              end
          end
        c = reshape(c, SIZE);
      end
      
      % ========================
    case {'struct' , 'parent'}
      c = feval(thisfunc, content, att_switch, struct([]), level+1);
      
      if ~prod(SIZE)==0
        c = reshape(c, SIZE);
      end
      
      if isfield(c, 'item') && strcmp(TYPE, 'struct')
        c = {c.item};          
      end
      
      % ========================
    case 'cell'
      tmp_c = feval(thisfunc, content, att_switch, {}, level+1);
      
      if ~prod(SIZE)==0
        tmp_c = reshape(tmp_c, SIZE);
      end
      
      if ~isempty(tmp_c)
        if isfield(tmp_c, 'item')
          c = {tmp_c.item};
        else
          % otherwise leave as is.
          c = tmp_c;
        end
      else
        c = {};
      end
      
      % ========================
    case 'sparse'
      S = feval(thisfunc, content, att_switch, {}, 1); % {} to ensure c is of type 'cell'
      if isempty(S)
        % we have an empty sparse matrix
        c = sparse(SIZE(1), SIZE(2));
      else
        c = sparse(S{1}, S{2}, S{3}, SIZE(1), SIZE(2));
      end
      
      % ========================
    case 'complex'
      C = feval(thisfunc, content, att_switch, [], 1);
      c = complex(C{1}, C{2});
      
      if ~prod(SIZE)==0
        c = reshape(c, SIZE);
      end
      
      % ========================  
    otherwise
      error([' Typeclass << ', TYPE, ' >> unknown to xml parser.']);
      
  end  
  
  % now c contains the content variable
  
  if isempty(X) && IDX==1 && level==0
    if strcmp(NAME, 'item')
      % s = '<item>aaa</item>'
      X = {};
      X(IDX) = {c};
    else
      % s = '<root>aaa</root>'
      X = c;
    end
    
  elseif isempty(X) && IDX==1 && level>0
    if strcmp(NAME, 'item')
      % s = '<root><item>bbb</item></root>'
      % s = '<root><item idx="1">a</item><item idx="2">b</item></root>'
      X = {};
      X(IDX) = {c};
    else
      % s = '<root><a>bbb</a></root>'
      X = setfield(X, {IDX}, NAME, c);
    end
    
  elseif isempty(X) && IDX>1 && level==0
    % s = '<root idx="4">hello</root>'
    % s = '<item idx="4">hello</item>'
    X = {}; 
    X(IDX) = {c};
    
  elseif isempty(X) && IDX>1 && level>0
    % s = '<root><ch idx="4">aaaa</ch></root>'
    % s = '<item><ch idx="4">aaaa</ch></item>'
    if strcmp(NAME, 'item')
      X = {};
      X(IDX) = {c};
    else
      X = setfield(X, {IDX}, NAME, c);
    end
    
  elseif ~isempty(X) && IDX==1 && level==0
    % s = '<item idx="3">aaa</item><item idx="1">bbb</item>'
    if strcmp(NAME, 'item')
      X(IDX) = {c};
    else
      if ~(nargin<3)
        % Example: a.b = 111; d = xml_parse(str, '', a);
        % this only works if both are structs and X is not empty
        if isempty(X) || ~(isa(X, 'struct') && isa(c, 'struct'))
          X = c;
        else
          % transfer all fields from c to X
          N = fieldnames(c);
          for n=1:length(N)
            X = setfield(X, {IDX}, N{n}, c.(N{n}));
          end
        end
      else
        % s = '<root idx="3">aaa</root><root idx="1">bbb</root>'
        % s = '<root>aaa</root><root>bbb</root>'
        % s = '<a><b>444</b></a><a><b>555</b></a>'
        errstr = sprintf(['XML string cannot have two ''root'' entries at root level! \n',...
            'Possible solution: Use ''item'' tags instead.']);
        error(errstr);
      end
    end
    
  elseif ~isempty(X) && IDX==1 && level>0

    if strcmp(NAME, 'item')
      % s = '<root><item idx="2">bbb</item><item idx="1">ccc</item></root>'
      X(IDX) = {c};
    else
      % s = '<root><a idx="2">bbb</a><a idx="1">ccc</a></root>'
      X = setfield(X, {IDX}, NAME, c);
    end
    % BUT:
    % s = '<root><a idx="2"><b>ccc</b></a><a idx="1">ccc</a></root>'
    % fails because struct a has different content!
    
  elseif ~isempty(X) && IDX>1 && level==0
  
    % s = '<item idx="1">a</item><item idx="2">b</item>'
    % s = '<item idx="1">a</item><item idx="2">b</item><item idx="3">c</item>'
    if isa(X,'char') 
      % s = '<item idx="1">a</item><item idx="2">b</item>'
      X = {X}; 
      %else (if not char) we would have eg the third entry as X
      %s = '<item idx="1">a</item><item idx="2">b</item><item idx="3">c</item>'
      %and do not need to take action
    end
    X(IDX) = {c};
    
  elseif ~isempty(X) && IDX>1 && level>0
  
    % s = '<root><item idx="1">a</item><item idx="2">b</item><item idx="3">c</item></root>'
    if strcmp(NAME, 'item')
      if isa(X,'char') 
        % s = '<root><item idx="1">a</item><item idx="2">b</item></root>'
        X = {X}; 
      end
      X(IDX) = {c};
    else
      % s = '<root><a>bbb</a><a>ccc</a></root>'
      X = setfield(X, {IDX}, NAME, c);
    end
    
  else
  
    disp('This case cannot be processed:')
    disp(['isempty(X) = ', num2str(isempty(X))])
    disp(['class(X)   = ', class(X)])
    disp(['class(c)   = ', class(c)])
    disp(['IDX        = ', num2str(IDX)])
    disp(['LEVEL      = ', num2str(level)])
    disp('Please contact the author m.molinari@soton.ac.uk!');
  end
  
  clear c;  
  i = i+1;
  
end

if nargout > 1
  xmlV = xmlVersion;
end

if level == 0
  % before we finally leave we should clean up variable
  % xmlTBVersion used as persistent in this function
  clear xmlTBVersion
end

return


% ==========================================================
function V = charunsubst(V)
% re-substitutes functional characters
% e.g. from '&lt;' to '<'.

f = strfind(V, '&amp;');  % &&&&&&&&&&
for i=1:length(f)
  V = [V(1:f(i)-1), '&', V(f(i)+5:end)];
  f = f-4;
end
f = strfind(V, '&lt;');  % <<<<<<<<<<
for i=1:length(f)
  V = [V(1:f(i)-1), '<', V(f(i)+4:end)];
  f = f-3;
end
f = strfind(V, '&gt;');  % >>>>>>>>>>
for i=1:length(f)
  V = [V(1:f(i)-1), '>', V(f(i)+4:end)];
  f = f-3;
end
f = strfind(V, '&apos;');  % ''''''''''
for i=1:length(f)
  V = [V(1:f(i)-1), '''', V(f(i)+6:end)];
  f = f-5;
end
f = strfind(V, '&quot;');  % " " " " " "
for i=1:length(f)
  V = [V(1:f(i)-1), '"', V(f(i)+6:end)];
  f = f-5;
end

% ==========================================================
function f = strfind2(longstr, str)
% find positions of occurences of string str in longstr
if size(longstr,2) < size(str,2)
  f=[];
  return
else
  f = strfind(longstr,str);
end
