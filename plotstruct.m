function [returntxt] = plotstruct (structree,varargin)
%PLOTSTRUCT plots recursively an overview of a structure in tree format and lists also size and some data fields
%
% function [returntxt] = plotstruct (structree)
% 
% Plots recursive an overview of a struct in tree format 
% and lists also size and some data fields.
% (Q&D implementation, not the fastest, but may help to get 
% an overview)
%
%                               (c) 2002, B. Kaeferstein
%                                 berthold@kaefersten.de

global childstruc; % save memory for recursive call

% handle recursive call
if nargin ==1
    lvl = 0;
    returntxt = [];
    datatree.ROOT_is = structree;
    structree = [];
    structree = datatree;
    datatree = [];
else 
    lvl = varargin{1};
    returntxt = varargin{2};
end

if ~isstruct(structree)
    return
end

% make the tree
s = '';
for cntr = [1:lvl]
    s = [s,char(32),char(32),char(32)]; % empty spaces to indent at levels
end    

fnames = fieldnames (structree);
if ~isempty (fnames) 
    for cntr = [1:size(fnames,1)]
%         childstruc = getfield(structree,fnames{cntr});
        tmp = [structree.(fnames{cntr})];
        if isstruct(tmp)
            childstruc = tmp(1);
        else
            childstruc = tmp;
        end
        % TABs for formatting -> makes it easy to copy result to Word etc.!
        outtxt = [s,'.',fnames{cntr}, char(9),class(childstruc),char(9),...
                num2str(smartsize(childstruc)),char(9),smartcontents(childstruc)]; 
        %disp(outtxt);
        returntxt = strvcat(returntxt,outtxt);
        returntxt = plotstruct (childstruc,lvl+1,returntxt);
    end
end    

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function smartsize = smartsize(obj);
% Get the size as vector if different from 1,1,1,1 etc.
sizevec = size (obj);
idx = find (sizevec ~=1);
if isempty(idx)
    smartsize = 1;
else
    smartsize = sizevec;
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function smartcontents = smartcontents(obj)
% Get some contents of the current object.
% Maybe I'll implement some search function later...
if  isa(obj,'char')
    smartcontents = obj (1:min(end,30));
elseif isa(obj,'double')
    smartcontents = obj (1:min(end,3));
elseif isa(obj,'single')
    smartcontents = obj (1:min(end,3));
elseif isa(obj,'numeric')
    smartcontents = obj (1:min(end,3));
elseif isa(obj,'int8')
    smartcontents = double(obj (1:min(end,3)));  
elseif isa(obj,'int16')
    smartcontents = double(obj (1:min(end,3)));  
elseif isa(obj,'int32')
    smartcontents = double(obj (1:min(end,3)));  
elseif isa(obj,'uint8')
    smartcontents = double(obj (1:min(end,3)));  
elseif isa(obj,'uint16')
    smartcontents = double(obj (1:min(end,3)));  
elseif isa(obj,'uint32')
    smartcontents = double(obj (1:min(end,3)));  
else
    smartcontents = [];
end

smartcontents = smartcontents(:)';
smartcontents = num2str(smartcontents);





