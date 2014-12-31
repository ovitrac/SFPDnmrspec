function y=uncell(x,depth,dim,removeempty)
%UNCELL remove cell imbrication until arbitrary depth long dimension dim
%   syntax: y=uncell(x [,depth,dim,removeempty])
%           x: cell
%       depth: maximum recursive depth (default = Inf)
%         dim: dimension (default = 1)
% removeempty: true to remove empty entry (default = false)
%
%   Example1: uncell({ {{{1}}} {{{'a'}}} {{{1:3}}} })
%             returns {1 'a' [1 2 3]}'
%   Example2: uncell({ {{{1}}} {{{'a' 'b'}}} {{{1:3}}} })
%             returns {{1} {'a' 'b'}' {[1 2 3]}}'

%MS 2.1 - 06/03/11 - INRA\Olivier Vitrac - rev. 06/07/11

%Revision history
%21/05/11 fix mixed row and column cell vectors
%22/05/11 fix empty x
%06/07/11 add removeempty

% arg check
if ~iscell(x), error('x must be cell'); end
if nargin<2, depth = []; end
if nargin<3, dim = []; end
if nargin<4, removeempty = []; end
if isempty(x), y=x; return; end
if isempty(depth), depth = Inf; end
if isempty(dim), dim = 1; end
if isempty(removeempty), removeempty = false; end
if depth==0 || ~all(cellfun(@(z) iscell(z),x(:))), y =x; return; end
siz = size(x{1});
if ~all(cellfun(@(z) all(size(z)==siz),x))
    if removeempty
        y = uncell(x(~cellfun('isempty',x)));
    else
        y=x; return
    end
end

% uncell
y = uncell(cat(dim,x{:}),depth-1,dim);