function [s,figposition] = StructDlgPrepare(pvalue,pname,pbounds,pflag,varargin)
%STRUCTDLGPREPARE convert structures into standard inputs for StructDlg UI controls 
% Final UI control styles and options in StructDlg are set with a cell array of up to four elements: 
% { value/definition labels limits protected_flag }
%
% SYNTAX 1: s = StructDlgPrepare(pvalue [,pname,pbounds,pflag])
%           s = StructDlgPrepare(pvalue [,pname,pbounds,pflag,'property1',value1,'property2',value2,...])
%   pvalue,pname,pbounds,pflag are structures with the same fields (only the fields matchine those of pvalue are considered)
%       pvalue.propertyname  = property value
%         pname.propertyname = real name (default = '')
%        pbounds.propertyname = [minvalue maxvalue] (default =[-Inf Inf])
%         pflag.propertyname = protection flag (default = false)
%
% SYNTAX 2: s = StructDlgPrepare(object)
%           s = StructDlgPrepare(object,[],[],[],'property1',value1,'property2',value2,...)
%   object = structure with fields
%   propertyname = structure('value',property value,...
%                            'name','real name',...
%                          'bounds',[minvalue maxvalue],...
%                            'flag', true or false)
% OPTION for SYNTAX 2: field length sets the size of the editable area in characters (only for text and numbers)
%
% For files (uigetfile/uiputfile) and directories (uigetdir), the meaning of 'value', 'name', 'bounds' and 'flag' are different
%   value = default file/directory name
%   name  = real name of the field
%  bounds = '*.ext1', or '*.ext1; *.ext2', etc. to set the possible extensions
%         = {'*.ext1'} or {'*.ext1; *.ext2'} etc. as above
%         = {'*.ext1' 'explanation for ext1'} or {'*.ext1' 'explanation for ext1' ; '*.ext2' 'explanation for ext2'} etc
%           to set an intelligible text besides the extension
%    flag = one of the following keyword 'getsinglefile', 'getmultiplefile', 'putsinglefile',  'getdir'
%        to indicate the kinf of file/directory selector to open
%       Note: for directories, bounds is not used
%
% List of implemented properties/values
%   'xmargin', 'ymargin' to set the margins along x and y (default value = 10)
%   'horizontalalignment' with possible values 'left', 'center' (default), 'right'
%   'verticalalignment' with possible values 'bottom', 'middle' (default), 'top'
%
%
% OUTPUTS
% s = StructDlgPrepare(...)
% [s,figposition] = StructDlgPrepare(...)
%   s = structure with fields
%       s.propertyname = {pvalue.propertyname pname.propertyname pbounds.propertyname pflag.propertyname}
%   figposition to be used with StructDlg
%
%
% SIMPLE EXAMPLE WITH SYNTAX 1
%   p = struct('number',12,'text','my text')
%   n = struct('number','this is a *','text','this a *')
%   b = struct('number',[0 20]);
%   l = struct('number',80);
%   [s,figposition]=StructDlgPrepare(p,n,b,[],l)
%   r=StructDlg(s,'My form to be filled out',[],figposition)
%
% SIMPLE EXAMPLE WITH SYNTAX 2 (object-oriented one)
%   p0 = struct(...
%       'number',struct('value',12,'name','this is *','bounds',[0 20],'flag',false),...
%       'text',struct('value','my text','name','this is a *','bounds',[],'flag',false) ...
%       );
%   [s0,figposition]=StructDlgPrepare(p0)
%   r0=StructDlg(s0,'My form to be filled out',[],figposition)
%
% SIMPLE EXAMPLE as above with a controlled length (only with SYNTAX 2)
%   p = struct(...
%       'number',struct('value',12,'name','this is *','bounds',[0 20],'flag',false),...
%       'text',struct('value','my text','name','this is a *','bounds',[],'flag',false,'length',50) ...
%       );
%   [s,figposition]=StructDlgPrepare(p,[],[],[],'horizontalalignment','right','verticalalignment','top')
%   r=StructDlg(s,'My form to be filled out',[],figposition)
% 
% AS ABOVE with simple file interface (only the first entry is used)
%  p.getsinglefile = struct('value','myxmlfile.xml','name','Choose a file to load','bounds',{{'*.xml','XML files'}},'flag','getsinglefile')
%  p.multiplefile = struct('value','','name','Pick several files','bounds',{{'*.xml','XML files'}},'flag','getmultiplefile')
%  p.putsinglefile = struct('value','myxmlfile.xml','name','Give a name to the file to save','bounds',{{'*.xml','XML files'}},'flag','putsinglefile')
%  p.getdir = struct('value',pwd,'name','Choose a directory','flag','getdir')
% [s,figposition]=StructDlgPrepare(p);
% r = StructDlg(s,'My form to be filled out',[],figposition)
%
% Generalized example with nested structure (SYNTAX 2 only)
% ATTENTION: do not forget the field value
% DO NOT USE UIGETFILE or UIGETDIR in child structures
%   p2= p;
%   p2.nested_struct = struct('value',p0);
%   p2.nested_struct.value.text.length = 120;
%   p2.nested_struct.value.text.value = 'this is my nested text';
%   p2.nested_struct.value.text.name = 'nested text';
%   [s2,figposition]=StructDlgPrepare(p2);
%   r = StructDlg(s2,'My form to be filled out',[],figposition)
%
% Generalized example with list
%   p1 = p2;
%   p1.list_of_elements = struct('value',{{'item1' 'item2' 'item3'}});
%   p1.nested_struct.value.nestedlist = struct('value',{{'item1' 'item2' 'item3' 'item4' 'item5'}},'name','my *','flag',true);  
%   [s,figposition]=StructDlgPrepare(p1)
%    r = StructDlg(s,'My form to be filled out',[],figposition)
%
%
% For debugging, close all open GUI figures
% deleteAllFigures = @(~, ~) delete(findobj(0,'type','figure')); h = figure; set(h,'DeleteFcn',deleteAllFigures);


% See also: StructDlg, argcheck
%
% Documentation of Structdlg
% web(fullfile(find_path_toolbox('migration'),'Structdlg/StructDlg.html'))

% Migration 2.2 -11/10/2014 - INRA\Olivier Vitrac - rev. 15/10/2014

% Revision history
% 14/10/2014 release candidate
% 15/10/2014 add figposition, implement lists

% Default values
strgetsinglefile = 'uigetfile({''%s'',''%s''},''%s'',''%s'')';
strgetmultiplefile = 'uigetfile({''%s'',''%s''},''%s'',''%s'',''MultiSelect'',''on'')';
strputsinglefile = 'uiputfile({''%s'',''%s''},''%s'',''%s'')';
strgetdir = 'uigetdir(''%s'',''%s'')';
default = struct('name','','limit',[-Inf +Inf],'flag',false,'length',[]);
getsinglefile = @(filename,pattern,patternname,title) {sprintf(strgetsinglefile,pattern,patternname,title,filename)};
getmultiplefile = @(filename,pattern,patternname,title) {sprintf(strgetmultiplefile,pattern,patternname,title,filename)};
putsinglefile = @(filename,pattern,patternname,title) {sprintf(strputsinglefile,pattern,patternname,title,filename)};
getdir = @(folder,title) {sprintf(strgetdir,folder,title)};
defaultpositions = struct(...
    'verticalalignment','middle',...
    'horizontalalignment','center',...
    'xmargin',10,...
    'ymargin',10 ...
);

% argcheck
if nargin<1, error('one argument is at least required'), end
if nargin<2, pname = struct([]); end
if nargin<3, pbounds = struct([]); end
if nargin<4, pflag = struct([]); end
if isempty(pvalue), s = struct([]); return, end
if ~isstruct(pvalue), error('the first argument (pvalue) must be a structure'), end
if all(structfun(@isstruct,pvalue)) && all(structfun(@(x) isfield(x,'value'),pvalue))
    for p=fieldnames(pvalue)'
        if ~isfield(pvalue.(p{1}),'name'), pvalue.(p{1}).name = ''; end
        if ~isfield(pvalue.(p{1}),'bounds'), pvalue.(p{1}).bounds = []; end
        if ~isfield(pvalue.(p{1}),'flag'), pvalue.(p{1}).flag = default.flag; end
        if ~isfield(pvalue.(p{1}),'length'), pvalue.(p{1}).length = default.length; end
    end
    plength = structfun(@(x) x.length,pvalue,'UniformOutput',false);
    pflag = structfun(@(x) x.flag,pvalue,'UniformOutput',false);
    pbounds = structfun(@(x) x.bounds,pvalue,'UniformOutput',false);
    pname = structfun(@(x) x.name,pvalue,'UniformOutput',false);
    pvalue = structfun(@(x) x.value,pvalue,'UniformOutput',false);
else
    plength = struct([]);
end
[windowposition,remainargs] = argcheck(varargin,defaultpositions);

if ~isempty(pname) && ~isstruct(pname), error('the second argument (pname) must be a structure'), end
if ~isempty(pbounds) && ~isstruct(pbounds), error('the third argument (pbounds) must be a structure'), end
if ~isempty(pflag) && ~isstruct(pflag), error('the fourth argument (pflag) must be a structure'), end
if ~isempty(plength) && ~isstruct(plength), error('the optional argument (plength) must be a structure'), end

% default values for all properties
p = fieldnames(pvalue); np = length(p);
pname_default = cell2struct(repmat({default.name},np,1),p);
pbounds_default = cell2struct(repmat({default.limit},np,1),p);
pflag_default = cell2struct(repmat({default.flag},np,1),p);
plength_default = cell2struct(repmat({default.length},np,1),p);
plength = argcheck(remainargs,plength);

% set default arguments
pname = argcheck(pname,pname_default);
pbounds = argcheck(pbounds,pbounds_default);
pflag = argcheck(pflag,pflag_default);
plength = argcheck(plength,plength_default);
islengthdefined = all(structfun(@isnumeric,plength)) && any(structfun(@length,plength));

% assembling the output
s = pvalue;
maxlength = 0;
maxlengthtext = 0;
for eachp = p'
    if isstruct(pvalue.(eachp{1}))
        s.(eachp{1}) = StructDlgPrepare(pvalue.(eachp{1})); %,pname.(eachp{1}),pbounds.(eachp{1}),pflag.(eachp{1}),plength.(eachp{1}));
    elseif ischar(pflag.(eachp{1}))
        filenameui = pvalue.(eachp{1});
        titleui    = pname.(eachp{1}); if isempty(titleui), titleui = ''; end
        if iscell(pbounds.(eachp{1}))
            patternui  = pbounds.(eachp{1}){1};
            if length(pbounds.(eachp{1}))>1, patternnfoui = pbounds.(eachp{1}){2}; else  patternnfoui = ''; end
        else 
            patternnfoui = pbounds.(eachp{1});
        end
        if isempty(patternnfoui), patternnfoui = sprintf('%s files',upper(patternui)); end
        switch pflag.(eachp{1})
            case 'getsinglefile'
                s.(eachp{1}) = {getsinglefile(filenameui,patternui,patternnfoui,titleui)};
                maxlength = max(maxlength,32);
            case 'getmultiplefile'
                s.(eachp{1}) = {getmultiplefile(filenameui,patternui,patternnfoui,titleui)};
                maxlength = max(maxlength,32);
            case 'putsinglefile'
                s.(eachp{1}) = {putsinglefile(filenameui,patternui,patternnfoui,titleui)}; 
                maxlength = max(maxlength,32);
            case 'getdir'
                s.(eachp{1}) = {getdir(filenameui,titleui)};      
                maxlength = max(maxlength,20);
            otherwise
                error('unknwown flag type ''%s'' in ''%s''',pflag.(eachp{1}),eachp{1})
        end
    elseif iscell(pvalue.(eachp{1}))
        if islengthdefined
            s.(eachp{1}) = {pvalue.(eachp{1}) pname.(eachp{1}) [] pflag.(eachp{1}) plength.(eachp{1})};
        else
            s.(eachp{1}) = {pvalue.(eachp{1}) pname.(eachp{1}) [] pflag.(eachp{1}) };
        end
    else
        ischartype = ischar(pvalue.(eachp{1})) && isempty(strfind(pvalue.(eachp{1}),'this.'));
        s.(eachp{1}) = { pvalue.(eachp{1}), pname.(eachp{1}) []  pflag.(eachp{1})};
        pbounds.(eachp{1}) = pbounds.(eachp{1})(~isnan(pbounds.(eachp{1})));
        if ~ischartype  && ~isempty(pbounds.(eachp{1}))
            if length(pbounds.(eachp{1}))==1, pbounds.(eachp{1}) = [0 pbounds.(eachp{1})]; end
            s.(eachp{1}){3} = pbounds.(eachp{1})(1:2);
        end
        if ischartype
            maxlengthtext = max(maxlengthtext,length(pvalue.(eachp{1})));
        end
    end
    if islengthdefined && iscell(s.(eachp{1}))
        s.(eachp{1}){5} = plength.(eachp{1});
        if ~isempty(plength.(eachp{1})), maxlength = max(maxlength,plength.(eachp{1})); end
    end
    
end

% guess best figposition based on the behavior StructDlg
if nargout>1
    [WidthScreen,HeightScreen] = get_screen_size('char');
    Height = np * 1.5 + 5;
    if maxlength==0 % || ~islengthdefined
        if maxlengthtext>0
            maxlength = maxlengthtext + 10;
        else
            maxlength = 20;
        end
    end
    Width = maxlength + 40;
    switch windowposition.horizontalalignment
        case 'left'
            xbox = windowposition.xmargin;
        case 'right'
            xbox = WidthScreen - Width - windowposition.xmargin;
        case 'center'
            xbox = round((WidthScreen-Width)/2);
        otherwise
            error('unknown horizontalalignment value = ''%s''',windowposition.horizontalalignment)
    end
    switch windowposition.verticalalignment
        case 'bottom'
            ybox = windowposition.ymargin;
        case 'top'
            ybox = HeightScreen - Height - windowposition.ymargin;
        case 'middle'
            ybox = round((HeightScreen-Height)/2);
        otherwise
            error('unknown verticalalignment value = ''%s''',windowposition.verticalalignment)
    end
    
    figposition = [xbox ybox Width Height];
end


%% PRINVATE FUNCTIONS
% GET_SCREEN_SIZE
function [WidthScreen,HeightScreen] = get_screen_size(units)
if ~nargin, units = 'pixels'; end
old_units = get(0,'units'); set(0,'units','pixels');
try
    monitorpos = get(0,'MonitorPositions');
    ismultiplemonitors = size(monitorpos,1)>1;
    if isunix
        ifirstmonitor = find(monitorpos(:,1)==0 | monitorpos(:,2)==0);
    else
        ifirstmonitor = find(monitorpos(:,1)==1 & monitorpos(:,2)==1);
    end
catch
    ismultiplemonitors = false;
end
set(0,'units',units);
if ismultiplemonitors && any(ifirstmonitor)
    monitorpos = get(0,'MonitorPositions');
    siz = monitorpos(ifirstmonitor,:);
else
    siz = get(0,'ScreenSize');
end
set(0,'units',old_units);
WidthScreen = siz(3);
HeightScreen = siz(4);