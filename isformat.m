function format_match = isformat(X,format,case_sense)
% ISFORMAT compare une chaîne à un ou une liste de format(s)
%		ex. format_match = isformat(X,format,[case_sense])
%		  format = string or cells of string
%           OLD SYNTAX for format: any expression with the following metacharacteurs
%               ?  = any character
%               #  = any digit
%               #* = any number (list of digits)
%               ~  = not a digit
%               @  = not a word character ~([a-z_A-Z0-9])              
%               *  = any list of characters
%               ==> any combinations and occurences are possible
%           NEW SYNTAX (including the keyword __REGEXP__) for format based on REGEXP
%               format = {'__REGEXP__',expr1,expr2,...}
%               see regexp for details
%           NB: without the keyword __REGEXP__, both syntaxes can be mixed (control the result)
%		  case_sense = 0 (by defauult)

% See: http://jspwiki.org/wiki/RegularExpressionSyntax
% See: web('jar:file:///C:/MATLAB73/help/techdoc/help.jar!/matlab_prog/f0-56452.html', '-helpbrowser')

% INRA\TCP-IP - 26/02/01 - Olivier Vitrac - rev. 04/12/05

% Revision History
% 04/12/05 use regexp instead of comp (more powerfull)
% 11/08/07 protect the symbol '.' in the regular expression (with '\.')
% 22/10/07 fix ?



% definition
case_sense_default = false;
keyword = '__REGEXP__';

% arg check
if nargin<3, case_sense = []; end
if isempty(case_sense), case_sense = case_sense_default; end

if ~iscell(format)
   format_match = isformat(X,{format},case_sense); %comp(X,format,case_sense);
else
    nformat = length(format);
    format_match = zeros(1,0);
    if strcmp(format{1},keyword) && nformat>1 % % new syntax        
        if case_sense
            res = regexp(X,format(2:end));
        else
            res = regexpi(X,format(2:end));
        end
        for i = 1:nformat-1
            if ~isempty(res{i})
                format_match(end+1) = i;
            end
        end
    else % old syntax
        if case_sense
            for i = 1:nformat
                if regexp(X,convertjockers(format{i})) % comp(X,format{i},case_sense);
                    format_match(end+1) = i; %format_match(end+1) = [format_match i];
                end
            end
        else
            for i = 1:nformat
                if regexpi(X,convertjockers(format{i})) % comp(X,format{i},case_sense);
                    format_match(end+1) = i; %format_match(end+1) = [format_match i];
                end
            end            
        end
    end
end

function xs = convertjockers(x)
% gateway function for capatibility with regexp
xs =  ['\<' x '\>'];
xs = strrep(xs,'.','\.');       % protect '.'
xs = strrep(xs,'?','.');        % any character
xs = strrep(xs,'#*','\d+');     % any number
xs = strrep(xs,'*','[^]*');     % any list of characters
xs = strrep(xs,'#','\d');       % any digit
xs = strrep(xs,'~','\D');       % not a digit
xs = strrep(xs,'@','\W');       % not a word character


% [OLD FUNCTIONS, NOT REQUIRED ANYMORE IN MATLAB 7.x]
% function res = comp(X,fmt,case_sense)
% lX = length(X);
% lf = length(fmt);
% indf = jocker_delete(fmt);
% indX = intersect(1:lX,indf);
% if case_sense
%    res = strcmp(X(indX),fmt(indf));
% else
%    res = strcmp(upper(X(indX)),upper(fmt(indf)));
% end
% u = find(fmt == '*');
% if any(u)
%    u = u(1);
%    if u>indf(1)
%    	i_u = find((u-indf == min(u-indf(find((u-indf>0)))))); % premier indice valide pointant avant le *
%    else, i_u = 0; end
%    if lX>=indf(end)
%         iX_succ = []; if_succ = iX_succ;
% 			%if i_u>0
%         %	if indf(i_u) < length(indf)
% 			%		indX = union(indX,indf(i_u+1):lX);
%         %   	if_succ = indf(i_u+1):indf(end);
%         %   	iX_succ = indX(end-length(if_succ)+1:end);
%         %   end
%         %else
%         if i_u<length(indf)
%            indX = union(indX,indf(i_u+1):lX);
%            if_succ = indf(i_u+1):indf(end);
%            iX_succ = indX(end-length(if_succ)+1:end);
% 			end
%         if case_sense
%            prec = strcmp(X(indf(1:i_u)),fmt(indf(1:i_u))); % validation de la chaîne gauche
%            succ = strcmp(X(iX_succ),fmt(if_succ));
%         else
%            prec = strcmp(upper(X(indf(1:i_u))),upper(fmt(indf(1:i_u)))); % validation de la chaîne gauche         
%            succ = strcmp(upper(X(iX_succ)),upper(fmt(if_succ)));
%         end
%       res = res | (prec & succ);
%    end
% else
%    res = res & (length(X)==length(fmt));
% end
%    
%    
% function ind = jocker_delete(X)
% ind = find(	(X ~= '#') & ... 
%    				(X ~= '~') & ...
% 					(X ~= '?') & ...
% 					(X ~= '*') & ...
% 					(X ~= '@'));