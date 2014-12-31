function [seqc,seqt] = grpseq(seq,step,formatnum,dash,comma,colon)
%GRPSEQ finds and groups consecutive numbers into cells
%   syntax   s = grpconsecutivenum(num,[step,format])
%        [s,t] = grpconsecutivenum(...)
%
%   Inputs:
%    num: mx1 or 1xm vector of integers (positive)
%   step: positive integer (default=1)
% format: number format (default='%d'), see sprintf
%         format of bounds can be set by using a 2x1 or 1x2 string array (see last examples for details)
%          3x1 or 1x3 can be used when step>1
%   dash: group operator (default='-')
%  comma: list separator (default=',')
%  colon: colon operator (default=':')
%
%   Outputs
%      s: nx1 cell array where s{i} is vector of consecutive integers
%      t: formatted string matching the sequence s (',' is the separator and '-' the consecutive operator)
%         Note: when step>1 '-' is replaced by ':step:'
%
%   Examples
%         [~,t]=grpseq([1:3 5 6:10])
%         t =
%         1-3,5-10
%
%
%         [s,t]=grpseq([1:2 5 6:10]*10,10), s{:}, t
%         s = 
%             [2x1 double]
%             [6x1 double]
%         ans =
%             10
%             20
%         ans =
%             50
%             60
%             70
%             80
%             90
%            100
%         t =
%         10,20,50:10:100
%
%
%         [~,t]=grpseq([1:3 5 6:10],[],'M%0.2d')
%         t =
%         M01-M03,M05-M10
%
%
%         [~,t]=grpseq([1:3 5 6:10],[],{'M%0.2d' '%0.2d'})
%         t =
%         M01-03,M05-10
%
%
%         [~,t]=grpseq([1:3 5 6:10],[],{'M_{%0.2d}' '_{%0.2d}'},'_-')
%         t =
%         M_{01}_-_{03},M_{05}_-_{10}
%
%         figure, title(t,'fontsize',16)

% MS 2.1 - 17/09/13 - INRA\Olivier Vitrac - rev. 23/10/2013

% revision history
% 23/10/2013 add format

% default
step_default = 1;
formatnum_default = '%d';
comma_default = ',';
dash_default = '-';
colon_default = ':';

% arg check
if nargin<1, error('one argument is required'), end
if nargin<2, step = []; end
if nargin<3, formatnum = ''; end
if nargin<4, dash = ''; end
if nargin<5, comma = ''; end
if nargin<6, colon = ''; end
if isempty(step), step = step_default; end
if isempty(formatnum), formatnum = formatnum_default; end
if isempty(dash), dash = dash_default; end
if isempty(comma), comma = comma_default; end
if isempty(colon), colon = colon_default; end
if ~isnumeric(seq), error('seq must be numeric'), end
seq = unique(round(seq(:)));
nseq = length(seq);
step = round(step(1));
if step == 1, formatmaxlength = 2; else formatmaxlength = 3; end %maximum size of accepted formatnum
if ischar(formatnum), formatnum = repmat({formatnum},1,formatmaxlength); end
if ~iscellstr(formatnum) || numel(formatnum)~=formatmaxlength
    error('format must be a char or a %dx1 o 1x%d cell array of strings',formatmaxlength)
end
formatsingle = ['%s%s' formatnum{1}];
formatinterval = ['%s%s' formatnum{1} '%s' formatnum{2}];
if step==1, formatcolon = ''; else formatcolon = [colon formatnum{3} colon]; end

% do the job
if ~isempty(seq)
    splitsize=diff(find(diff([-inf;seq-(1:nseq)'*step;inf])));
    seqc=mat2cell(seq,splitsize,1);
else
    seqc = {};
end

% text
if nargout>1
    nc = length(seqc);
    seqt = '';
    for i=1:nc
        if i>1, sep1 = ','; else sep1=''; end
        if step==1
            if length(seqc{i})==2, sep2=comma; else sep2=dash; end
        else
            if length(seqc{i})==2, sep2=comma; else sep2=sprintf(formatcolon,step); end
        end
        if length(seqc{i})>1
            seqt = sprintf(formatinterval,seqt,sep1,seqc{i}(1),sep2,seqc{i}(end));
        else
            seqt = sprintf(formatsingle,seqt,sep1,seqc{i});
        end
    end 
end