function [hout,propout] = plotpub(x,y,varargin)
%PLOTPUB extend 1D PLOT capabilities for scientific publishing (initially based on PLOTEASY)
%  By default, PLOTPUB(X,Y,'-o') plots continuous curves with sparse symbols.
%  PLOTPUB(X,Y,'-a') works also with arbitrary symbols not included in PLOT.
%  PLOTPUB(X,Y,'-pa') introduces an extended color scheme with shortcuts for brown, pink, gray
%  PLOTPUB(X,{Y1 Y2 Y3},{'-o' '-a'},'markersize',[8 12]) accepts any combination of data and properties (with circular permutation)
%  PLOTPUB(X,Y,P) accepts multiple plots and multiple properties with simplified syntax and reusable properties P.
%  PLOTPUB(X,Y,[P1 P2]) combines several properties
%  PLOTPUB(...,'disp','pretty') improves the readibility of sparse curves
%  PLOTPUB(...,'type','sparse') is the default type for curves
%  PLOTPUB(...,'xlim',[xmin xmax]) modify properties to the current axes(implemented properties include xlim/ylim, xscale/yscale, xtick/ytick...)
%  >> PLOTPUB does not obey to figure properties: 'DefaultAxesColorOrder' and 'DefaultAxesLineStyleOrder'.
%
%   PLOTTING syntaxes
%  ====================
%        as PLOT        : plotpub(x,y,[s,'property1',value1,property2,value2...])
%                         plotpub(y,[s,'property1',value1,property2,value2...])
%      generalizing PLOT: plotpub(x,{y1 y2...},[s,'propertynum',[num1 num2...] ,'property2',{'string1' 'string2'...},...])
%                         plotpub({x1 x2...},{y1 y2...},[s,'propertynum',[num1 num2...] ,'property2',{'string1' 'string2'...},...])
%     reusing properties: plotpub(x,y,prop,[...])
% handles and properties: [h,prop] = plotpub(...);
%
%   PROPERTIES syntaxes
%  ====================
%     default properties: prop = plotpub;
%    creating properties: prop = plotpub('propertynum',[num1 num2...] ,'property2',{'string1' 'string2'...},...]):
%   mofifying properties: prop = plotpub(prop,'propertynum',[num1 num2...] ,'property2',{'string1' 'string2'...},...]):
%concatenation of properties: prop = plotpub([prop1 prop2],[...]):
% distribution of properties: propdist = plotpub(prop,[...]), where propdist is a nyx1 structure array
%
%   INPUTS
%       x,y can be cell arrays or matrices
%           Each element of the cell array is display with personalized styles.
%           A matrix is displayed with a same style
%           A matrix can be converted (here columnwise) into a cell array using num2cell(x,1).
%          >> 
%       x{i} or x: mx1 or mxn arrays 
%       y{i} or y: mxn arrays (use num2cell(y,1) to break a mxn y matrix into a 1xn cell array of mx1 vectors)
%       s: character string or string arrays
%          >> accept all PLOT syntaxes
%          >> non PLOT commands are interpreted as symbols (can be any ASCII characters)
%          >> since some characters can be missinterpreted, the following conventions are applied:
%             ->priority to color
%               'o' stands for 'orange' (priority to color')
%               'oo' stands for an orange symbol 'o'
%             ->non conventional symbols as non PLOT interpretable syntaxes and are plotted with TEXT)
%       prop: style properties coded as a structure, with a general format
%             prop.('property') = {'string1' 'string2'...} for any non numeric or matrix properties
%             prop.('property') = [value1 value2...] for any numeric properties
%             non numeric properties are: 'color' 'disp' 'fontname' 'interp' 'linestyle' 'marker' 'markeredgecolor'
%             'markerfacecolor' 'plotter' 'sign' 'tex' 'textcolor' 'textlinestyle' 'xscale'
%             numeric properties are: 'linewidth' 'markersize' 'nmarker' 'textlinewidth'
%          >> see the section RECOGNIZED PROPERTIES for details
%          >> missing properties are replaced by default values (defaultvalues = plotpub;)
%          >> non recognized properties generate an error
%          >> concatenation of properties such as [prop1 prop2] is accepted as inputs
%          >> prop = plotpub([prot1 prot2]) generates a property structure where all fields are concatenated
%          >> propdist = plotpub(prot) distributes in a structure array the properties of several plots,
%             where protdist(1), protdist(2)... are the properties for the different curves.
%       [...] parameter (lower case)/value pairs to specify additional properties (cell are authorized)
%             By convention, parameter values overlay always shortcuts.
%
%   OUTPUTS
%       h: vector of structure, which contains the handles of plotted objects
%           h(i).leg = vector of handles corresponding to xi and yi, which are comptabible with LEGEND (i.e. related to hidden objects)
%           h(i).line = corresponding vector of handles of possible line objects (i.e. visible)
%           h(i).marker = corresponding vector of handles of possible marker objects (i.e. visible)
%           h(i).text = corresponding vector of handles of possible text objects (i.e. visible)
%       All objects are tagged, with 'plotpub:type', where type = 'hidden', 'line', 'marker', 'text'
%         >> to legend specific plots, use: plotpub(h(subindex),...)
%         >> to legend different plots, concatenate the corresponding handes, example: plotpub([h1 h2],...)
%         >> plots from different axes or figures can be legended with legendpub (not possible with LEGEND)
%       prop: structure with fields matching the properties (see prop above)
%         >> prop.('property') is a 1xn vector or cell array, where n is numel(y)
%         >> when y is missing (PROPERTIES syntax), n is the largest length of the provided properties
%
%   RECOGNIZED PROPERTIES (Note: they are all vectorized, only lowercase names are accepted)
%       Properties compatible with plot
%           'color':      'k' or gray(3) or {'k' [0 1 1]} or any combination
%           'linestyle':  '-' or {'-' '--'} or any cell of strings or '' to remove a line
%           'linewidth':  [1 2] or any vector
%           'marker':     '*' or {'*' '^'} or any cell of acceptable symbols (see plot) or '' to remove a symbol
%                         ascii characters are also acceptable (they are plotted with TEXT with fontsize defined by MARKERSIZE)
%           'markersize': [4 8] or any vector
%           'markerfacecolor': as 'color' (note that 'none' stands for no color and nor for [1 1 1])
%           'markeredgecolor': as 'color
%       Text properties for non PLOT symbol properties (are replaced by markers ones when possible):
%           'tex' any tex string or combination of strings (default = '\rm'), an empty or 'none' string inactivate the TEXT interpreter
%           'fontname': 'helvetica' (default), 'times', 'Courier', 'FixedWidth'...
%           'markersize' replaces 'fontsize'
%           'markerfacecolor' replaces 'backgroundcolor' (default = 'none')
%           'markeredgecolor' replaces 'edgecolor' (default = 'none')
%           'textcolor' text color (default = 'k')
%                       >> implemented using the tex command \color[rgb]{R G B}
%           'textlinestyle'   set the linestyle of the rectangle's edge line (default = '-')
%           'textlinewidth'   set the linewith of the rectangle's edge line (default = 0.3)
%           'horizontalalignment' default = 'center'
%           'verticalalignment' default = 'middle'
%       Other properties (which are not standard in MATLAB)
%           'nmarker': number of requested markers (default = 12), a combination of values is possible
%           'nmarkermin': number of minimum marker (valid with disp='optimized...'), default = 3;
%           'xlimmarker': valid range for marker (default = [-Inf +Inf]);
%           'markerpriority': set how conflicts are resolved when disp='... ends'
%                       'bulk' (default): bulk points are preserved
%                       'ends': ends are preserved
%                       other values, do nothing
%           'plotter': plotter to be used (possibly a user or ad-hoc 1D plotter)
%                      'plot': Matlab plotter (default)
%                      'plotequi': as plot but remove lines between not equispaced points (no interpolation)
%                      'plotsmoothbar': plots bars as a normal 1D x,y plot
%                      'stairs': Matlab plotter
%                      'csaps': plot(c,fnval(csaps(x,y),x))
%                      'csapszz': plot(c,fnval(csaps(x,y,zz/100),x))
%                      >> a combination of plotter is possible {'plot' 'plotequi'}
%           'sign'   : 'all' (default) '+' (only positive values) '-', a combination is possible {'all' '+' '+' '-'}
%           'interp' : 'index' 'nearest','linear' (default),'cubic' or any combination
%           'disp'   : 'pretty' (default) prevents the superposition of symbols when a same xscale is used for all x{i}
%                      'pretty with ends' as 'pretty' but plots also the ends
%                      'optimized' precalculate the best common scale to display symbols
%                      'optimized with ends' as 'optimized' but plots also the ends
%                      'none'   no optimization
%           'type'  :  'sparse' (default), 'continuous'
%                      if 'sparse', 'o-' is interpreted as a continuous line with sparse symbols
%                      otherwise 'o-' is interpreted as a continuous line with continuous symbols
%       Axes properties (empty = default behavior)
%           'xlim'   : [xmin xmax] (use +/-Inf or for automatic value)
%                  >> default = not applied
%                  >> precedence rule: xmin = min of all xmin, max = max of all xmin
%           'ylim'   : [ymin ymax] idem
%           'xscale' : 'linear' (default) or 'log' or 'sqrt', a combination is possible {'linear' 'log' 'log'}
%                  >> sparse plots are updated according to xscale
%                  >> a 'log' value have a precedance over any 'linear' value
%           'yscale' : idem
%           'box'    : 'on' (default) or 'off'
%                  >> 'on' have a precedance over any 'off'
%           'xtick'  : vector of xtick values
%                  >> default = not applied
%                  >> precedence rule: combine all values
%           'ytick'  : vector of ytick values
%           'xticklabel': 'auto' (default), 'none' (removed the label)
%           'yticklabel': 'idem
%       Extented color scheme (MATLAB+SIMULINK+new personalized ones)
%           r: {'red'  [1 0 0]}
%           g: {'green'  [0 1 0]}
%           y: {'yellow'  [1 1 0]}
%           m: {'magenta'  [1 0 1]}
%           b: {'blue'  [0 0 1]}
%           k: {'black'  [0 0 0]}
%           w: {'white'  [1 1 1]}
%           G: {'gray'  [0.5000 0.5000 0.5000]}
%           c: {'cyan'  [0 1 1]}
%           f: {'darkGreen'  [0 0.5000 0]}
%           o: {'orange'  [1 0.7500 0]}
%           B: {'lightBlue'  [0.6800 0.9200 0]}
%           W: {'brown'  [0.6000 0.2000 0]}
%           p: {'pink'  [1 0.6000 0.7800]}
%           v: {'violet'  [0.5000 0.2000 0]}
%
%    See also: LEGENDPUB, DELETEPLOTPUB, GETPROP, SETPROP, SUBPLOTS, OBJ2IM, PLOTOBJ, PRINT_JPG, REFSLOPE, ADDZPLOTPUB
%
% EXAMPLES to show how to replace advantegeously the matlab PLOT function with the following general syntax:
%  plot(x1,y1,'string 1',x1,y2,'string 2','property1',value1,'property2',value2)
%    x1 = linspace(0,10,1000)';
%    x2 = linspace(0,10,1500)';
%    y1 = ncx2pdf(x1,4,2);
%    y2 = chi2pdf(x2,4);
%    figure, plotpub(x1,y1,'o-')
%    disp('it looks like plot but it displays an orange line without symbols ''o''')
%    figure, plotpub(x1,y1,'oo-')
%    disp('plots a continuous line with 12 symbols ''o''')
%    figure, plotpub(x1,y1,'ko-')
%    disp('plots a black line with 12 symbols ''o''')
%    figure, plotpub({x1 x2},{y1 y2},{'y-s' 'B:p'},'linewidth',[2 4],'markersize',[10 12],'markerfacecolor',{'c' [0 1 0]})
%    disp('a more sophisticated example')
%    figure, plotpub(y1,'marker','o','linestyle',':')
%    disp('rapid syntax')
%    figure, plotpub(magic(3),'interp','linear')
%    disp('note that the default ''interp'' value, ''cubic'' does not provide satisfactory value')
%    figure, plotpub(x1,{y1 y1/2},{'r-o' 'gs-'})
%    disp('a rapid combination of several properties')
%    figure, plotpub({1:3 1:4},{magic(3) magic(4)},{'ra-' 'cb-'},'interp','linear')
%    disp('an impossible figure with plot')
% 
% A LAST EXAMPLE
%  [h,p] = plotpub({1:3 1:3 1:6},{(1:3) 2*(1:3) 0.5*(1:6)},...
%      'color',{'W' 'b'},'marker',{'a' 'b' 'c'},'markeredgecolor','none',...
%       'markersize',18,'nmarker',[3 5],'textcolor',{'y' 'f'})
% To reuse the same properties for similar plots, use:
%   h2 = plotpub({1:3 1:3 1:6},{(1:3) 2*(1:3) 0.5*(1:6)},p);
% To display a legend
% legendpub(h2,{'leg1',{'leg1_line1' 'leg1_line2' 'leg1_line3' 'leg1_line4' 'line5'}},[],[],'interpreter','none')

% MS 2.1 - 09/08/07 - Olivier Vitrac - rev. 08/03/11

% history
% 09/08/07 Release candidate
% 12/08/07 major fixes, s introduced, help is rewritten
% 14/08/07 fix xscale=='log', add 'disp' property, additional help on x and y
% 17/08/07 fix color shortcuts, add properties syntax definition, additional help
% 18/08/07 fix 'pretty' for a mixture such as 'k-o' with 'ko' (only 'k-o' lines require to be optimized, fix hold
% 18/08/07 fix interpolation for matrices with repeated x values
% 19/08/07 add axes properties and distribution of properties into a structure array
% 21/08/07 fix prop when y is missing (allow prop as single input)
% 27/08/07 fix 'nearest'when NaN are generated, fix distribution of prop when nprop>2, add 'plotsmoothbar' as plotter
% 29/08/07 x and y can be 2D cell arrays (all LENGTH(x or y) are replaced by NUMEL(x or y))
% 30/08/07 fix plotpub([0 0],[0 1]) (it generated an errors since sparce cannot be interpolated)
% 30/08/07 xtick,ytick,xlim and ylim are not anymore parsed when they are provided as array (expected behavior)
% 05/09/07 same fix as above for counting ny when only properties are provided
% 08/09/07 add 'stairs' as plotter
% 09/09/07 add 'horizontalalignment', 'verticalalignment', 'xlimmarker', 'nmakermin' + fix empty(xi)
% 10/09/07 fix parsing of ex: plotpub('linewidth',[.3 .6 .3]) => 3 props and not 1 (required for refslope)
% 07/09/08 add 'csaps'and 'csasps' as plotter
% 02/12/08 add 'sqrt' as xscale
% 20/03/10 fix ny for colors
% 08/03/11 fix text symbols on a sqrt scale, accept mixed colors rgb and letter codes

% Default (configuration section, must be modified with care since parsing rules relies on the content of default)
default =    struct( ... parsing rules are indicated as comments
    'color',{{'k'}},... cell required
    'linestyle',{{'-'}},... cell required
    'linewidth',1,... numeric required
    'marker',{cellstr(fliplr('*<>vx+sohpd^.')')'},... cell required
    'markersize',8,... numeric required
    'nmarker', 12,... numeric required
    'nmarkermin',3,... numeric required
    'markerpriority',{{'bulk'}},... cell required 
    'xlimmarker',{{[]}},... no default value
    'markerfacecolor',{{'none'}},... cell required (and so on)
    'markeredgecolor',{{'k'}},...
    'textcolor',{{'k'}},...
    'textlinestyle',{{'-'}},...
    'textlinewidth',.2,...
    'tex',{{'\rm'}},...
    'fontname',{{'helvetica'}},...
    'horizontalalignment',{{'center'}},...
    'verticalalignment',{{'middle'}},...
    'plotter',{{'plot'}},...
    'sign',{{'all'}},...
    'xscale',{{'linear'}},...
    'interp',{{'cubic'}},...
    'disp',{{'pretty'}},...
    'type',{{'sparse'}},...
    'box',{{'on'}},...
    'yscale',{{''}},... no default value. Note {[]} as content means no parsing
    'xlim',{{[]}},...   no default value. Note {''} does not have the behavior as {[]} since '' is not numeric
    'ylim',{{[]}},...   no default value
    'xtick',{{[]}},...  no default value
    'ytick',{{[]}},...  no default value
    'xticklabel',{{''}},...  no default value
    'yticklabel',{{''}} ...  no default value
                     );
proplist = fieldnames(default)';
markerlist = '.ox+*sdv^<>ph';
tagid = 'plotpub:';
defaultmarker = 'x'; % for legend when a character is used
color2fullname = struct(...
    'r',{{'red' [1 0 0]}},...
    'g',{{'green' [0 1 0]}},...
    'y',{{'yellow', [1 1 0]}},...
    'm',{{'magenta' [1 0 1]}},...
    'b',{{'blue' [0 0 1]}},...
    'k',{{'black' [0 0 0]}},...
    'w',{{'white' [1 1 1]}},...
    'G',{{'gray' [0.5 0.5 0.5]}},...
    'c',{{'cyan' [0 1 1]}},...
    'f',{{'darkGreen' [0 0.5 0]}},...
    'o',{{'orange' [1 0.75 0]}},...
    'B',{{'lightBlue' [0.68 0.92 0]}},...
    'W',{{'brown' [0.6 0.2 0]}},...
    'p',{{'pink' [1 0.6 0.78]}},...
    'v',{{'violet' [0.5 0.2 0]}}...
    );
validlinestyle = {'*:*' '*-.*' '*--*' '*-*'}; % recognized linestyle (the order is important)
validcolor = fieldnames(color2fullname)'; % only the color letter codes are used (full color names are not accepted in this version)
xismissing = false;   % flag, true if x is omitted (as in plot)
yismissing = false;   % flag, true if y is omitted (use to retrieve properties objects)
plotshortcut = false; % flag, true if syntax similtat to PLOT(X,Y,S) (see plot for details)
effnargin = nargin; % since nargin is a function (not a variable) and it cannot be modified
providedprop = cell(1,0); % provided properties

% INPUT PARSING (almost compatible with PLOT syntax)
% no input to parse
if effnargin<1, hout = default; return, end % returns default parameters
% distribution of properties as a structure array
if effnargin==1 && isstruct(x) && numel(x)==1 %length(x)==1
    prop = x;
    proplist = fieldnames(prop)'; nproplist = length(proplist);
    nprop=Inf; for f=proplist; nprop = min(nprop,length(prop.(f{1}))); end
    proptmp = cell(nprop,nproplist);
    for i=1:nproplist, proptmp(:,i) = mat2cell(prop.(proplist{i})(1:nprop),1,ones(1,nprop)); end
    hout = cell2struct(proptmp,proplist,2);
    return
end
% parse x and y values
if ischar(x) || isstruct(x) || (iscell(x) && ischar(x{1}))   % no data at all (no x and no y)
    yismissing = true;
    if effnargin>1, varargin(2:end+1) = varargin; varargin{1}=y;end
    varargin(2:end+1) = varargin; varargin{1}=x; effnargin=effnargin+2;
else % parse at least y
    if effnargin==1, y = x; xismissing = true; effnargin=effnargin+1; end
    if effnargin>=2
        if ischar(y) || isstruct(y) || (iscell(y) && ischar(y{1})) %%elseif iscell(y) && ischar(y{1}), xismissing = true; end
            xismissing = true;
        end
        if xismissing && nargin>1
            varargin(2:end+1) = varargin; varargin{1}=y;
            y=x; effnargin=effnargin+1; % change the number of input arguments
        else
            prop = [];
        end
    end
    if xismissing % add missing x
        if isnumeric(y)
            if size(y,1)>1, x = (1:size(y,1))'; else x = 1:size(y,2); end
        elseif iscell(y)
            my = numel(y); %length(y);
            x = cell(1,my);
            for i=1:my
                if isnumeric(y{i})
                    if size(y{i},1)>1, x{i} = (1:size(y{i},1))'; else x{i} = 1:size(y{i},2); end
                else error('unrecognized y data, y{%d} is expected to be numeric or of type cell',i')
                end
            end
        else error('unrocognized y data, y is expected to be numeric or a cell')
        end
    end
end
% parse plot properties (shortcurts, structure or pair properties/values)
if effnargin>=3
    if isstruct(varargin{1}) % if it is a property object
        prop = varargin{1};
        if length(prop)>1         % check the size of prop
            proptmp=struct;
            for f=fieldnames(prop)'
                if ~ismember(f{1},proplist), error('the field ''%s'' does not exist',f{1}), end
                expectedtype = class(default.(f{1}));
                if ~isa(prop(1).(f{1}),expectedtype)
                    error('the field ''%s''is expected to be of class ''%s'',\nthe provided class is ''%s''',f{1},expectedtype,class(prop(1).(f{1})))
                end
                proptmp.(f{1}) = [prop.(f{1})]; % concatenation of all properties
            end
            prop = proptmp;
        end
    elseif ~ismember(varargin{1},proplist) % if it is not a recognized property
        prop = shortcuts(varargin{1},validcolor,validlinestyle,markerlist); % try to apply plot shortcuts (string syntax)
        if effnargin>3 && isstruct(varargin{2}) % shortcut and structure (redundant input)
            proptmp = varargin{2};
            for f=fieldnames(proptmp)'
                prop.(f{1}) = proptmp.(f{1}); % priority is given to additional parameters
            end
            varargin = varargin([1 3:end]); % remove the redundant input
            effnargin = effnargin -1;
        end % redundant
    else % no property
        prop = [];
    end
    if ~isempty(prop) && effnargin>3 % parse addinional pair parameters/values
        if mod(effnargin-3,2), error('parameters combined with shortcuts must be associated by pairs'), end
        proptmp = cell2struct(varargin(3:2:end),varargin(2:2:end-1),2);
        for f=fieldnames(proptmp)'
            prop.(f{1}) = proptmp.(f{1}); % priority is given to additional parameters
        end
    end
end
if isempty(prop) && mod(effnargin-2,2)>0, error('parameters must be given by pair parameter/value(s)'), end
if ~isstruct(prop) && effnargin>2, prop = cell2struct(varargin(2:2:end),varargin(1:2:end-1),2); end
if isempty(prop), prop = default; end
providedprop = intersect(fieldnames(prop),proplist); % list all provided properties


% check x and y values
if ~yismissing
    nx = numel(x); ny = numel(y); %nx = length(x); ny = length(y);
    if ~iscell(x), if size(x,1)==1, x=x'; end, nx = size(x,2); x = num2cell(x,1); end
    if ~iscell(y), if size(y,1)==1, y=y'; end, ny = size(y,2); y = num2cell(y,1); end
    if nx<ny, x = num2cell(repmat(x{1},1,ny),1); end
    if ny<nx, y = num2cell(repmat(y{1},1,nx),1); ny=nx; end
else % only ny is required
    ny = 0;
    for f=providedprop %fieldnames(prop)'
        if (isnumeric(prop.(f{1})) && iscell(default.(f{1})) && ~isempty(default.(f{1}){1}) ) || iscell(prop.(f{1}))
            if any(strfind(f{1},'color')), ny = max(ny,size(prop.(f{1}),1));
            else ny = max(ny,length(prop.(f{1}))); end
        elseif isnumeric(prop.(f{1})) && isnumeric(default.(f{1})) && length(default.(f{1}))==1
            ny = max(ny,length(prop.(f{1}))); % example 'linewidth'
        elseif ~isempty(prop.(f{1}))
            ny = max(ny,1);
        end
    end
end

% Check all properties
for f = fieldnames(prop)'; if ~isfield(default,f{1}), error('''%s'' is an invalid property',f{1}), end, end % check properties name
for f = proplist
    if ~isfield(prop,f{1}) || isempty(prop.(f{1})), prop.(f{1}) = default.(f{1}); end % add missing properties
    if any(strfind(f{1},'color')) % any color field
        if isnumeric(prop.(f{1}))
            if size(prop.(f{1}),2)~=3, error('rgb values of ''%s'' should be a nx3 matrix',f{1}), end
            if any((prop.(f{1})>255) | (prop.(f{1})<0)), error('rgb values of ''%s'' are ranged between 0 and 1 (or between 0 and 255)',f{1}), end
            if any((prop.(f{1})>1)), prop.(f{1}) = prop.(f{1})/255; end
            prop.(f{1}) = mat2cell(prop.(f{1}),ones(size(prop.(f{1}),1),1),3)'; 
        end
    end
    if isnumeric(default.(f{1})) && ~isnumeric(prop.(f{1})), error('the property ''%s'' must be numeric',f{1}), end
    if iscell(default.(f{1})) && ~iscell(prop.(f{1})) % check the appropriate type of each property
        if isnumeric(default.(f{1}){1})
            if isempty(default.(f{1}){1})
                prop.(f{1}) = {prop.(f{1})}; % no parsing (used for xlim,ylim,xtick,ytick...)
            else
                prop.(f{1}) = num2cell(prop.(f{1}),2); % parsing
            end
        else
            prop.(f{1}) = cellstr(prop.(f{1}));
        end
    end 
    prop.(f{1}) = prop.(f{1})(mod(0:ny-1,length(prop.(f{1})))+1);
end
if yismissing
    if nargout>1, error('only the properties are returned'), end
    hout = prop;
    return
end

% current axes properties
currentax = gca; currenthold = ishold;
if ~currenthold, delete(get(currentax,'children')); end
hold on

% check the number of lines
nylines = zeros(ny,1);
for i=1:ny
    if ~isempty(prop.linestyle{i}) && ~strcmp(prop.linestyle{i},'none') && ...
            (any(strfind(prop.disp{i},'pretty')) || any(strfind(prop.disp{i},'optimized'))) && ...
            strcmp(prop.type{i},'sparse')
        nylines(i:end) = nylines(end)+1;
    end
end

% global xscale optimization to display symbols
if length(unique(prop.nmarker))~=1
    prop.disp = repmat({'none'},1,ny); % pretty is enough
    disp(sprintf('''disp'' property is downgraded to ''none''\nsince different ''nmarker'' values are applied'))
end
if ismember('optimized',prop.disp) || ismember('optimized with ends',prop.disp)
    if (ny>1) && (nx>1)
        minxi = Inf; maxxi = -Inf;
        for i=1:nx
            minxi = min(minxi,min(x{i}(:)));
            maxxi = max(maxxi,max(x{i}(:)));
        end
        nxxi = max(prop.nmarker);
        if ismember('log',prop.xscale); % logscale
            if maxxi<0, error('all data are negative, no possible log scale for x'), end
            if minxi<0, error('some data are negative, the ''disp'' property cannot be set to ''optimized'''), end
            maxxi = log10(maxxi); minxi = log10(minxi);
            dxopt = (maxxi-minxi)/(nxxi+1);
            xopt = logspace(minxi,maxxi-dxopt,nxxi)';
            set(currentax,'xscale','log')
            forcedlogscale = true; % force log scale
        else % a linear scale is assumed
            dxopt = (maxxi-minxi)/(nxxi+1);
            xopt = linspace(minxi,maxxi-dxopt,nxxi)';
            forcedlogscale = false;
        end
    else
        prop.disp = repmat({'pretty'},1,ny); % pretty is enough
        disp('''disp'' property is downgraded to ''pretty''')
        disp('since a similar xscale is used')
    end
end

% Plots
h = repmat(struct('leg',[],'line',[],'marker',[],'text',[]),1,ny);
for i=1:ny % for all x{i} y{i}
    % check input sizes
    if size(x{i},1)==1, x{i}=x{i}'; end
    if size(y{i},1)==1, y{i}=y{i}'; end
    nxx = size(x{i},2); nyy = size(y{i},2);
    if nxx>1 && (nxx~=nyy), error('incompatible sizes for for x{%d} and y{%d}',i,i), end
    % sign properties are on the combination of all columns
    switch prop.sign{i} % keep only valid values
        case 'all', j = (1:size(y{i},1))';
        case '+',   j = find(prod(single(y{i}>0),2));
        case '-',   j = find(prod(single(y{i}<0),2));
    end
    if any(prop.linestyle{i}) && ~strcmp(prop.linestyle{i},'none') % a line is required
        if strcmpi(prop.plotter{i},'plot')
            h(i).line = plot(x{i}(j,:),y{i}(j,:),prop.linestyle{i},...
                                       'linewidth',prop.linewidth(i),...
                                       'color',rgbconverter(prop.color{i},color2fullname,'k'));
        elseif any(strfind(lower(prop.plotter{i}),'csaps')) %strcmpi(prop.plotter{i},'csaps')
            lambda = strread(lower(prop.plotter{i}),'csaps%d');
            if any(lambda)
                sp = csaps(x{i}(j,:),y{i}(j,:),lambda/100);
            else
                [sp,lambda] = csaps(x{i}(j,:),y{i}(j,:));
                dispf('lambda = %0.3g%%',lambda*100)
            end
            h(i).line = plot(x{i}(j,:),fnval(sp,x{i}(j,:)),prop.linestyle{i},...
                                       'linewidth',prop.linewidth(i),...
                                       'color',rgbconverter(prop.color{i},color2fullname,'k'));            
        elseif strcmpi(prop.plotter{i},'plotequi')
            h(i).line = plotequi(x{i},y{i},prop.linestyle{i},...
                'linewidth',prop.linewidth(i),...
                'color',rgbconverter(prop.color{i},color2fullname,'k'));
        elseif strcmpi(prop.plotter{i},'plotsmoothbar')
            h(i).line = plotsmoothbar(x{i},y{i},prop.linestyle{i},...
                'linewidth',prop.linewidth(i),...
                'color',rgbconverter(prop.color{i},color2fullname,'k'));
        elseif strcmpi(prop.plotter{i},'stairs')
            h(i).line = stairs(x{i}(j,:),y{i}(j,:),prop.linestyle{i},...
                                       'linewidth',prop.linewidth(i),...
                                       'color',rgbconverter(prop.color{i},color2fullname,'k'));            
        else
            error('error in plotting [x{d},y{d}]: ''%s'' is an invalid plotter.\n It must be implemented before use',i,i,prop.plotter{i})
        end
        set(h(i).line,'Tag',[tagid 'line'])
        h(i).leg = h(i).line;
        if any(prop.marker{i}) && ~strcmp(prop.marker{i},'none') && prop.markersize(i)>0 % markers are required
            % X SCALE
            ax = get(currentax,'xlim'); %xlim;
            notoptimized = true;
            if any(strfind(prop.disp{i},'optimized')) % optimized for all data
                if forcedlogscale
                    if strcmp(prop.disp{i},'optimized with ends')
                         xi = [xopt(1);xopt*10^(dxopt*(i-1)/ny+dxopt/2);xopt(end)*10^dxopt];
                         xi = endsconflictsolver(xi,dxopt,prop.markerpriority{i},'log');
                    else
                        xi = xopt*10^(dxopt*(i-1)/ny+dxopt/2);
                    end
                else
                    if strcmp(prop.disp{i},'optimized with ends')
                        xi = [xopt(1);xopt+(dxopt*(i-1)+dxopt/2)/ny;xopt(end)+dxopt];
                        xi = endsconflictsolver(xi,dxopt,prop.markerpriority{i},'linear');
                    else
                        xi = xopt+(dxopt*(i-1)+dxopt/2)/ny;
                    end
                end
                xi = xi( (xi>=min(x{i}(:))) & (xi<=max(x{i}(:))) );
                if length(xi)>=prop.nmarkermin(i)
                    notoptimized = false;
                else
                    disp(sprintf('(x{%d},y{%d}): x range too small (%d<%d) to be used with ''%s''\n switch to ''%s''',i,i,length(xi),prop.nmarkermin(i),prop.disp{i},strrep(prop.disp{i},'optimized','pretty')))
                end
            end
            if notoptimized % non globally optimized ('pretty' or ''pretty with ends')
                if strcmp(prop.xscale{i},'linear') || strcmp(prop.xscale{i},'sqrt')
                    minxi = min(x{i}(j));
                    maxxi = max(x{i}(j));
                    nxxi = round(prop.nmarker(i)*size(y{i},1)/length(j));
                    if (ny==1) || ~any(strfind(prop.disp{i},'pretty')), dx=0; else dx = (maxxi-minxi)/(nxxi+1); end
                    if strfind(prop.disp{i},'ends') %strcmp(prop.disp{i},'pretty with ends')
%                         xi = [minxi;linspace(minxi,maxxi-dx,nxxi)'+(dx*(i-1)+dx/2)/ny;maxxi];
                        xi = [minxi;linspace(minxi,maxxi-dx,nxxi)'+(dx*(nylines(i)-1)+dx/2)/max(1,nylines(end));maxxi];
                        xi = endsconflictsolver(xi,dx,prop.markerpriority{i},'linear');
                    else 
%                         xi = linspace(minxi,maxxi-dx,nxxi)'+(dx*(i-1)+dx/2)/ny;
                        xi = linspace(minxi,maxxi-dx,nxxi)'+(dx*(nylines(i)-1)+dx/2)/max(1,nylines(end));
                    end
                elseif strcmp(prop.xscale{i},'log')
                    u = find(x{i}(j)>0); % only positive values
                    j = j(u);
                    nxxi = round(prop.nmarker(i)*size(y{i},1)/length(j));
                    minxi = min(log10(x{i}(j)));
                    maxxi = max(log10(x{i}(j)));
                    if ny==1 || ~any(strfind(prop.disp{i},'pretty')), dx=0; else dx = (maxxi-minxi)/(nxxi+1); end
                    if strfind(prop.disp{i},'ends') %strcmp(prop.disp{i},'pretty with ends')
%                         xi = [10^minxi;logspace(minxi,maxxi-dx,nxxi)'*10^(dx*(i-1)/ny+dx/2);10^maxxi];
                        xi = [10^minxi;logspace(minxi,maxxi-dx,nxxi)'*10^(dx*(nylines(i)-1)/max(1,nylines(end))+dx/2);10^maxxi];
                        xi = endsconflictsolver(xi,dx,prop.markerpriority{i},'log');
                    else 
%                         xi = logspace(minxi,maxxi-dx,nxxi)'*10^(dx*(i-1)/ny+dx/2);
                        xi = logspace(minxi,maxxi-dx,nxxi)'*10^(dx*(nylines(i)-1)/max(1,nylines(end))+dx/2);
                    end
                    set(currentax,'xscale','log')
                else
                    error('unable to interpret the scale ''%s'' for y{%d}',prop.scale{i},i)
                end
            end
            % INTPERPOLATION (if sparse)
            if strcmp(prop.type{i},'sparse')
                [xtmp,yi] = interpinternal(x{i}(j,:),y{i}(j,:),xi,prop.interp{i});
                if any(xtmp)
                    xi = xtmp;
                    if nxx>1
                        if ~strcmp(prop.interp{i},'nearest')
                            warning('REDUNDANCY detected: only the first column of y{%d} is used.\nUse cell2num(y) as input to remove this warning.',i) %#ok<WNTAG>
                        else
                            warning('interpolation method =''interp'': only the first column of y{%d} is used.\nUse cell2num(y) as input to remove this warning.',i) %#ok<WNTAG>
                        end
                    end
                end
                % filtering according 'xlimmarker' (added 9/9/7)
                if length(prop.xlimmarker{i})==2
                    indfilter = (xi>=prop.xlimmarker{i}(1)) & (xi<=prop.xlimmarker{i}(2));
                    xi = xi(indfilter); yi = yi(indfilter);
                elseif ~isempty(prop.xlimmarker{i})
                    disp(sprintf('(x{%d},y{%d}): unable to interpret ''xlimmarker''',i,i))
                end
            else % no interpolation
                xi = x{i}(j,:);
                yi = y{i}(j,:);
            end   
                
%             % PLOTS
            if length(prop.marker{i})==1 && any(markerlist==prop.marker{i}) %  PLOT marker
                h(i).leg = plot(xi,yi,[prop.marker{i} prop.linestyle{i}],...
                    'linewidth',prop.linewidth(i),...
                    'color',rgbconverter(prop.color{i},color2fullname),...
                    'markersize',prop.markersize(i),...
                    'markeredgecolor',rgbconverter(prop.markeredgecolor{i},color2fullname),...
                    'markerfacecolor',rgbconverter(prop.markerfacecolor{i},color2fullname) ...
                    );
                set(h(i).leg,'visible','off','Tag',[tagid 'leg'])
                h(i).marker = plot(xi,yi,prop.marker{i},...
                    'markersize',prop.markersize(i),...
                    'markeredgecolor',rgbconverter(prop.markeredgecolor{i},color2fullname),...
                    'markerfacecolor',rgbconverter(prop.markerfacecolor{i},color2fullname) ...
                    );
                set(h(i).marker,'Tag',[tagid 'marker'])
            else % TEXT marker
                h(i).leg = plot(xi,yi,[prop.linestyle{i}],'linewidth',prop.linewidth(i),'color',rgbconverter(prop.color{i},color2fullname,'k'));
                set(h(i).leg,'visible','off','Tag',[tagid 'leg'])
                if ~isempty(prop.tex{i}) && ~strcmp(prop.tex{i},'none')
                    interpreter = 'tex';
                    texstring = [prop.tex{i} texcolor(rgbconverter(prop.textcolor{i},color2fullname))];
                else
                    interpreter = 'none'; texstring = '';
                end
                for iyy=1:nyy
                    h(i).text(:,iyy) = text(xi,yi(:,iyy),[texstring prop.marker{i}],...
                        'fontsize',prop.markersize(i),...
                        'linestyle',prop.textlinestyle{i},...
                        'linewidth',prop.textlinewidth(i),...
                        'edgecolor',rgbconverter(prop.markeredgecolor{i},color2fullname),...
                        'backgroundcolor',rgbconverter(prop.markerfacecolor{i},color2fullname),...
                        'horizontalalignment',prop.horizontalalignment{i},...
                        'verticalalignment',prop.verticalalignment{i},...
                        'interpreter',interpreter);
                end
                set(h(i).text,'Tag',[tagid 'text'])
            end % plot marker
        end
    else % not a line
        if (length(prop.marker{i})==1) && any(markerlist==prop.marker{i}) % a valid marker
            h(i).marker = plot(x{i}(j),y{i}(j),'marker',prop.marker{i},...
                'markeredgecolor',rgbconverter(prop.markeredgecolor{i},color2fullname),...
                'markerfacecolor',rgbconverter(prop.markerfacecolor{i},color2fullname),...
                'markersize',prop.markersize(i),...
                'linestyle',prop.linestyle{i});
            h(i).leg = h(i).marker;
            set(h(i).marker,'Tag',[tagid 'marker'])
        elseif ~strcmp(prop.marker{i},'none') && ~isempty(prop.marker{i})  % a text
            h(i).leg = plot(x{i}(j),y{i},defaultmarker); axis tight, ax = axis;
            set(h(i).leg,'visible','off','Tag',[tagid 'leg'])
            if ~isempty(prop.tex{i}) && ~strcmp(prop.tex{i},'none')
                interpreter = 'tex';
                texstring = [prop.tex{i} texcolor(rgbconverter(prop.textcolor{i},color2fullname))];
            else
                interpreter = 'none'; texstring = '';
            end
            ixx=0;
             for iyy=1:nyy
                 ixx = min(ixx+1,nxx);
                 h(i).text(:,iyy) = text(x{i}(j,ixx),y{i}(j,iyy),[texstring prop.marker{i}],...
                     'fontsize',prop.markersize(i),...
                     'edgecolor',rgbconverter(prop.markeredgecolor{i},color2fullname),...
                     'backgroundcolor',rgbconverter(prop.markerfacecolor{i},color2fullname),...
                     'linestyle',prop.textlinestyle{i},...
                     'linewidth',prop.textlinewidth(i),...
                     'horizontalalignment',prop.horizontalalignment{i},...
                     'verticalalignment',prop.verticalalignment{i},...
                     'interpreter',interpreter);
             end
            set(h(i).text,'Tag',[tagid 'text'])
            set(currentax,'xlim',ax(1:2),'ylim',ax(3:4))
        end % a valid marker
    end % a line
end % each i
set(gca,'box','on')

% Axes properties
if ismember('on',prop.box), set(currentax,'box','on'), else set(currentax,'box','off'), end
if ismember('log',prop.xscale), set(currentax,'xscale','log'), elseif ismember('linear',prop.xscale), set(currentax,'xscale','linear'), end
if ismember('log',prop.yscale), set(currentax,'yscale','log'), elseif ismember('linear',prop.yscale), set(currentax,'yscale','linear'), end
ax  = struct('xlim',[NaN NaN],'ylim',[NaN NaN],'xtick',[],'ytick',[]);
for i=1:ny % apply precedence
    for f={'xtick' 'ytick' 'xlim' 'ylim'}
        if ~isnumeric(prop.(f{1}){i}), error('the property ''%s'' related to y{%d} must numeric or []',f{1},i), end
    end
    for f={'xlim' 'ylim'}
        if any(prop.(f{1}){i})
            if length(prop.(f{1}){i})~=2, error('the property ''%s'' related to y{%d} must be [xmin xmax]',f{1},i), end
            ax.(f{1})(1) = min(min(prop.(f{1}){i}),ax.(f{1})(1));
            ax.(f{1})(2) = max(max(prop.(f{1}){i}),ax.(f{1})(2));
        end
    end
    for f={'xtick' 'ytick'}
        if any(prop.(f{1}){i}), ax.(f{1}) = union(ax.(f{1}),prop.(f{1}){i}); end
    end
end
if any(ax.xlim) && any(ax.xtick)
    ax.xlim(1) = min(ax.xlim(1),min(ax.xtick));
    ax.xlim(2) = max(ax.xlim(2),max(ax.xtick));
end
if any(ax.ylim) && any(ax.ytick)
    ax.ylim(1) = min(ax.ylim(1),min(ax.ytick));
    ax.ylim(2) = max(ax.ylim(2),max(ax.ytick));
end
for f={'xtick' 'ytick' 'xlim' 'ylim'}
    if any(ax.(f{1}))
        tst = diff(ax.(f{1}));
        if any(tst) && ~any(isnan(tst)) && ~all(tst==0)
            set(currentax,f{1},ax.(f{1}))
        end
    end
end
if ismember('none',prop.xticklabel), set(currentax,'xticklabel',' '); end
if ismember('none',prop.yticklabel), set(currentax,'yticklabel',' '); end

% sqrt xscale
if ismember('sqrt',prop.xscale)
    xtick_sqrt      = sqrt(get(currentax,'xtick'));
    xticklabel_sqrt = get(currentax,'xticklabel');
    xlim_sqrt       = sqrt(get(currentax,'xlim'));
    for i = 1:ny
        for f=fieldnames(h(i))'
            if strcmp(f{1},'text') && ~isempty(h(i).text)
                for eachtext = 1:length(h(i).text)
                    xtext = get(h(i).text(eachtext),'position');
                    xtext(1) = sqrt(xtext(1));
                    set(h(i).text(eachtext),'position',xtext)
                end
            elseif ~strcmp(f{1},'leg')
                set(h(i).(f{1}),'Xdata',sqrt(get(h(i).(f{1}),'Xdata')))
            end
        end
    end
    set(currentax,'xtick',xtick_sqrt,'xticklabel',xticklabel_sqrt,'xlim',xlim_sqrt)
end


% Outputs
if ~currenthold, hold off, end
if nargout>0, hout = h; end
if nargout>1, propout = prop; end
drawnow


% ===================
%  private functions
% ===================
function s = shortcuts(str,colors,linestyles,markers)
% Parser for shortcuts
if ~ischar(str) || iscell(str), s = []; end
if ~iscell(str), str = {str}; end, nstr = length(str);
isashortcut = true;
s = struct('linestyle',{{}},'marker',{{}},'color',{{}},'markeredgecolor',{{}},'markerfacecolor',{{}},'textcolor',{{}},'textlinestyle',{{}}); % prefetch
for i=1:nstr % scan all shortcuts (cell array)
    if ~ischar(str{i}), error('plot properties must be a string or a structure'), end
    str{i} = str{i}(:)';
    % search for a valid linestyle (STEP 1)
    j = isformat([' ' str{i} ' '],linestyles); % the last version of isformat is required (some inconsistency may occur with the old version)
    if any(j)
        linestyle = linestyles{(j(1))}(linestyles{(j(1))}~='*');  % always the first
        s.linestyle{i} = linestyle;
        k = strfind(str{i},linestyle);
        str{i} = [str{i}(1:k-1) str{i}(k+length(linestyle):end)]; % remove the current linestyle
    else
        s.linestyle{i} = 'none';
    end
    % scan all letter color code (STEP 2)
    ncolors = length(colors);
    foundcolors = ''; foundposition = [];
    for j=1:ncolors
        k = strfind(str{i},colors{j});
        if k, foundcolors(end+1) = colors{j}; foundposition(end+1)=k(1); end
    end
    if any(foundcolors) % if any color code
        [k,j] = min(foundposition); % take only the first one
                s.color{i} = foundcolors(j);
                str{i} = [str{i}(1:k(1)-1) str{i}(k(1)+1:end)]; % remove the current color code
        foundcolors = foundcolors(j);
    end
    % search for a marker (if it is a letter code => marker otherwise the remaining is a text) (STEP 3)
    nchar = length(str{i});
    foundmarkers = []; foundposition = [];
    if nchar>1
        istext = true; % no letter code possible, it is text
    elseif nchar==1
        k = find(str{i}==markers);
        if any(k)
            istext = false;
            s.marker{i} = markers(k);
            if foundcolors
                s.markeredgecolor{i} = foundcolors;
                s.markerfacecolor{i} = 'none';
            end
            istext = false;
        else
            istext = true;
        end
    else
        istext = false;
    end
    if istext
        s.marker{i} = str{i};
        s.markeredgecolor{i} = 'none';
        s.markerfacecolor{i} = 'none';
        s.textlinestyle{i} = 'none';
        if foundcolors, s.textcolor{i}=foundcolors; end     
    end
    if length(s.marker)<i || isempty(s.marker{i}), s.marker{i}='none'; end
end % next i

function rgb = rgbconverter(color,color2fullname,none)
% return the appropriate rgb value
if nargin<3, none = 'none'; end
if iscell(color), color = color{1}; end % after 8/3/11
if iscell(color), error('invalid color properties ''%s''',color{:}), end
% if iscell(char), error('invalid color properties ''%s''',color), end % before 8/3/11
if isempty(color), rgb = none; return, end
if ischar(color)
    if strcmp(color,'none')
        rgb = none;
    else
        rgb = color2fullname.(color(1)){2};
        if isempty(rgb), rgb = none; end
    end
else
    rgb = color(1,1:3);
end

function tex = texcolor(rgb)
% convert any rgb value into an appropriate tex command
if isempty(rgb) || ischar(rgb)
    rgb = '';
else
    tex = sprintf('\\color[rgb]{%0.3g %0.3g %0.3g}',rgb);
end

function [xout,yi] = interpinternal(x,y,xi,method)
% internal interpolation function, which accepts redundancy
% if redundancy detected or method = 'nearest' only the first column of y is considered
extrapval = NaN;
nx = size(x,2);
[my,ny] = size(y);
mxi = length(xi);
if (nx>1) || (ny>1)
    yi = zeros(mxi,ny);
    xout = [];
    i=1;
    while isempty(xout) && i<=ny
        [xout,yi(:,i)] = interpinternal(x(:,min(i,nx)),y(:,i),xi,method); % partial recursion for efficiency
        i = i+1;
    end
    if any(xout), yi = yi(:,1); end
else
    innan = find((~isnan(x))&(~isnan(y))); % all NaN are removed;
    x = x(innan);
    y = y(innan);
    dx = diff(x);
    if strcmp(method,'nearest') || any(dx==0) % redundancy
        [xunique,j] = unique(x); nxunique = length(xunique);
        if nxunique>1
            k = interp1(xunique,(1:nxunique)',xi,'nearest','extrap');
            k = k((~isnan(k)) & (k>=1) & (k<=nxunique)); % remove NaN and invalid extrap values
            xout = xunique(k);
            yi = y(j(k));
        else % no available x values for interpolation
            xout = xi;
            yi = interp1(1:length(y),y,linspace(1,length(y),length(xi)),method);
        end
    else % not redundant
        xout = [];
        yi = interp1(x,y,xi,method,extrapval);
    end
end

function x=endsconflictsolver(x,dx,method,scale)
% solve conflict when ends are too close (<dx/2) to the next point
if length(x)>2
    if strcmp(scale,'log'), close = log10([x(2)/x(1) x(end)/x(end-1)])<dx;
    else close = [x(2)-x(1) x(end)-x(end-1)]<dx; end
    if close(1)
        if strcmp(method,'bulk'), x = x(2:end);
        elseif strcmp(method,'ends'), x = x([1 3:end]); end
    end
    if close(2)
        if strcmp(method,'bulk'), x = x(1:end-1);
        elseif strcmp(method,'ends'), x = x([1:end-2 end]); end
    end
end