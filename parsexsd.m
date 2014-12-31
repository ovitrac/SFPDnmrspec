function X = parsexsd(filename,parsingrule)
%  PARSEXSD parse a XSD file
%   Syntax: molstruct = parsexsd(filename)
%   Options: molstruct = parsexsd(filename,'autoupdate') 
%
% X=parsexsd('C:\Data\Olivier\INRA\MS\TH 4\Migrants\Nonanal Files\Documents\nonanal.xsd')
% X=parsexsd('C:\Data\Olivier\INRA\MS\TH 4\Migrants_Database_eg_Piringer\Database_Piringer Files\Documents\P158 Disco Min\P158.xsd')
% X=parsexsd('C:\Data\Olivier\INRA\MS\TH 4\Polymères\PP Files\Documents\Polypropylene.xsd')
%
% See also: LOADXSD HIS2ARC ARCREAD

% MS-MATLAB 1.0 - 11/02/04 - Olivier Vitrac - rev. 18/04/12

% Revision History
%   01/03/04 fix the case of an empty xsd files
%   18/09/07 add id file
%   12/11/08 ass "see also"
%   18/04/12 fix fileparts
%   18/04/12 recompilation http://www.artefact.tk/software/matlab/xml/
%            cd(fullfile(find_path_toolbox('MS'),'@xmltree/private'))
%            mex -O xml_findstr.c

% Definitions
% -----------
% keyword at atomic scale {'xsd format' 'MS format'}
% first argument is assumed to be ID
atomkw = { 'Atom3d'   'ato'
           'Bond'    'bnd'   };
% properties at atomic scale {'xsd format' 'MS format' 'type' 'conversion'}
% first argument is assumed to be ID
% valid type are 'num' 'lst' and 'str'
atoprop{1} =     { 'ID'             'Ai'      'num'
                   'XYZ'            'Axyz'    'num'
                   'Connections'    'Ac'      'lst'
                   'Charge'         'AC'      'num'
                   'ForcefieldType' 'AF'      'str'
                   'Components'     'A'       'str'
                   'Name'           'id'      'str'
                 };
atoprop{2} =     { 
                    'ID'             'Bj'      'num'
                    'Connects'       'B'      'lst'
                 };
iatoprop = {[] []};
% properties at molecular scale {'xsd format' 'MS format' isnumeric}
% first argument is assumed to be ID
molprop = { 'ID'        'Mi'      'num'
            'XYZ'       'Mxyz'  'num'
            'name'      'M'     'str'
        };
molist = [];
numformatprop = '%f,'; % sscanf properties



% arg check
[pathstr,name,ext] = fileparts(filename);
if isempty(pathstr), pathstr = cd; end
if isempty(ext), ext = '.xsd'; end
fullfilename = fullfile(pathstr,[name ext]);
if ~exist(fullfilename,'file'), error(sprintf('the file ''%s'' does not exist'),fullfilename), end
if nargin<2, parsingrule = 'autoupdate'; end
autoupdate = strcmpi(parsingrule,'autoupdate');

% load xml file
fxsd = xmltree(fullfilename);
lxsd = length(fxsd);

% parser
firstmol = 1;
firstatom = ones(size(atoprop));
imolecule = 0;
iatom = [];
X = [];
for i=1:lxsd
    name = get(fxsd,i,'name');
    % check for Atom3D or Bond keys
    j = find(ismember(atomkw,name));
    if any(j) % if valid key
        % parent molecule id
        idmoli = parent(fxsd,i);                    % internal id
        moli = attributes(fxsd,'get',idmoli,1);     % MS id
        moli = str2num(moli.val);
        if ~ismember(moli,molist) % if new molecule found
            molist = [molist;moli]; % add molecule to list of molecules
            imolecule = length(molist);
            moliprop = attributes(fxsd,'get',idmoli); % list of prop
            moliprop = [moliprop{:}]; % vectorization
            if firstmol % if never parsed before (define rule for parsing)
                [found,pos] = isformat_multiple(molprop(:,1)',{moliprop.key});
                if any(~found)
                    for k=find(~found), disp(sprintf('the molecular property ''%s'' is missing',molprop{k,1})), end
                    error('... see above')
                end
                imolprop = [pos{:}];
                if ~autoupdate,  firstmol = 0; end
            end
            for k = 1:length(found); %size(molprop,1)
                switch molprop{k,3}
                    case 'num'
                        val = sscanf(moliprop(imolprop(k)).val,numformatprop)';
                        if isempty(val)
                            error(sprintf('missing value for the property ''%s''',molprop{k,1}));
                        end
                        X(imolecule).(molprop{k,2}) =  val;
                    otherwise
                        X(imolecule).(molprop{k,2}) =  moliprop(imolprop(k)).val;
                end
            end
            iatom(imolecule,1:length(atoprop)) = 0;
            iatom(imolecule,j) = 1;
        else % is a listed molecule
            imolecule = find(molist==moli);
            iatom(imolecule,j) = iatom(imolecule,j)+1;
        end % new molecule
        %                 % parse information at atomistic scale
        atomiprop = attributes(fxsd,'get',i);   % list of prop
        atomiprop = [atomiprop{:}];             % vectorization
        if firstatom(j) % define parsing rule
            [found,pos] = isformat_multiple(atoprop{j}(:,1)',{atomiprop.key});
            if length(found)<size(atoprop{j},1)
                notfound = setdiff(1:size(atoprop{j},1),found);
                for k=notfound
                    disp(sprintf('Warning: the atomic property ''%s'' is missing',atoprop{j}{k,1}))
                    atomiprop(end+1).val = 'NaN';
                    pos{end+1} = length(atomiprop);
                end
                % reordering
                [allitems,idx] = sort([found notfound]);
                pos = pos(idx); 
            end
            iatoprop{j} = [pos{:}];
            if ~autoupdate, firstatom(j) = 0; end
        end
        for k = 1:size(atoprop{j},1)
            switch atoprop{j}{k,3}
                case 'num'
                    val = sscanf(atomiprop(iatoprop{j}(k)).val,numformatprop)';
                    if isempty(val)
                        disp(sprintf('Warning: bad parsing for the property ''%s''',atoprop{j}{k,1}))
                    else
                        X(imolecule).(atoprop{j}{k,2})(iatom(imolecule,j),:) =  val;
                    end
                case 'lst'
                    val = sscanf(atomiprop(iatoprop{j}(k)).val,numformatprop)';
                    X(imolecule).(atoprop{j}{k,2}){iatom(imolecule,j),1} =  val;
                otherwise
                    X(imolecule).(atoprop{j}{k,2}){iatom(imolecule,j),1} = atomiprop(iatoprop{j}(k)).val;
            end
        end  
    end % any(j)
end % for i


