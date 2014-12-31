function combout = createcombination(entity,class,nchoose)
%CREATECOMBINATION generate all combinations when choosing n entities between m total entities (n<m) with 1 constraint (acceptable combinations have at least one entity belongging to each class)
% SYNTAXES  combout = creatrcombination(entity,class)
%           combout = creatrcombination(entity,class,5)
% INPUTS
%          ENTITY: m x 1 array of entities introduced
%           CLASS: m x 1 array of classes that belongging entities
%         NCHOOSE: numeric (nchoose < m), number of entities chose between 'm' entities (default)
% OUTPUTS
%         COMBOUT: k x nchoose matrix, 'k' combinations when choosing 'nchoose' entities
%
% EXAMPLE 1
% entities = (1:10)'; class = [1 1 2 3 4 5 6 6 6 7]';
% combout = createcombination(entities,class,3);
% combout =
% 
%      1     2     3
%      1     3     4
%      2     3     4
% EXAMPLE 2
% entities = (1:10)'; class = [1 2 2 3 4 5 6 6 6 7]';
% combout = createcombination(entities,class,2);
 
% RMNSPEC v 0.1 - 05/12/2013 - INRA\Olivier Vitrac, LNE\Mai Nguyen - rev. 29/12/2013

% Revision history
% 29/12/2013 fix nchoose, add recursion

% default
nchoose_default = 1;

% argcheck
if nargin < 1, error('At least one argument is required'), end
if nargin < 2, class = []; end
if nargin < 3, nchoose = []; end
if isempty(nchoose), nchoose = nchoose_default; end

% recursion
n_nchoose = length(nchoose);
if n_nchoose>1
    combout = cell(n_nchoose,1);
    for i=1:n_nchoose
        combout{i} = createcombination(entity,class,nchoose(i));
    end
    return
end

nentity = length(entity);
if isempty(class), class = 1:nentity; end
[~,orderentity] = sort(entity);
classofeachentity = class(orderentity);

% guess all possible combinations
comb = @(n) nchoosek(entity(1:find(class<=n,1,'last')),n);

% generate combinations and check that all classes are present in combinations 
combout = comb(nchoose);
ncomb = size(combout,1);
ok = false(ncomb,1);
requiredclass = 1:class(nchoose);
for i=1:ncomb
    ok(i) = all(ismember(requiredclass,classofeachentity(combout(i,:))));
end
combout = combout(ok,:);

