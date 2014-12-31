function txt=formatsci(num,varargin)
%FORMATSCI print numbers as scientific Tex strings
%   Syntax: txt=formatsci(num [,param1,value1,param2,value2])
%
%   Pair properties
%       power10: anonymous function to calculate the power of 10, default value = @(x) floor(log10(x))
%    texpattern: valid formatting tex string with sprintf placeholders, default value =  %0.3g\cdot10^{%d}
%                Note that '\' will be automatically protected
%       pattern: alternative pattern to be compared when economoy is used
%   Keywords
%       txt = formatsci(...,'economy')
%       
%   Example: formatsci([0.5e-15 12 1 1000 100000,-45,-89000,1e2],'economy') gives:
%               '5\cdot10^{-16}'    '12'    '1'    '10^{3}'    '10^{5}'    '-45'    '-8.9\cdot10^{4}'    '100'
%            without 'economy' the result is:
%               '5\cdot10^{-16}'    '1.2\cdot10^{1}'    '10^{0}'    '10^{3}'    '10^{5}'    '-4.5\cdot10^{1}'    '-8.9\cdot10^{4}'    '10^{2}'


% MS 2.1 - 25/10/11 - INRA\Olivier Vitrac - rev. 08/12/13


% revision history
% 09/01/12 fix 0, add negative numbers
% 27/01/12 add pattern and economy
% 29/03/12 remove 10^0 fro results
% 08/12/13 improve output

% default values
param_default = struct(...
    'power10',@(x) floor(log10(x)),...
    'texpattern','%0.3g\cdot10^{%d}',...
    'pattern','%0.4g' ...
    );
keywordlist = 'economy';

% arg checks
param = argcheck(varargin,param_default,keywordlist);
param.texpattern = regexprep(param.texpattern,'\\','\\\');
texwhenremeq1 = uncell(regexp(param.texpattern,'^.*(10.*)','tokens'));

% power 10 and remainder
sgn = sign(num);
num = sgn.*num;
p10 = param.power10(num);
rem = num./10.^(p10);

% print collected values
% txt = sprintf(param.texpattern,[rem;p10]);
% txt = arrayfun(@(x,y) sprintf(param.texpattern,x,y),rem,p10,'UniformOutput',false);
n = numel(num);
txt = cell(1,n);
for i=1:n
    if num(i)==0
        txt{i} = '0';
    elseif abs(rem(i)-1)<sqrt(eps)/1e3;
        txt{i} = sprintf(texwhenremeq1{1},p10(i));
    elseif abs(p10(i))~=0
        txt{i} = sprintf(param.texpattern,rem(i),p10(i));
    else
        txt{i} = sprintf(param.pattern,rem(i));
    end
    if param.economy
        tmp  = sprintf(param.pattern,num(i));
        tmptex = regexprep(txt{i},{'10\^','\^|{|}|\\cdot'},{'e',''});
        if (length(tmp)<=length(tmptex)) && ~any(tmp=='e'), txt{i} = tmp; end
    end
    if sgn(i)<0, txt{i}=['-' txt{i}]; end
end
if n<2, txt = txt{1}; end