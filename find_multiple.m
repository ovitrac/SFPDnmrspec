function [ind,val_trouve,val_non_trouve] = find_multiple(rech,base,mise_en_forme)
% FIND_MULTIPLE recherche les occurences des éléments d'un sous-vecteur dans un autre vecteur
%		ex. ind = find_multiple(rech,base)
%		options : [ind,val_trouve,val_non_trouve] = find_multiple(rech,base,mise_en_forme)
%
%		mise_en_forme = 0 (défaut), 1 ou 2
%
%			-> [i,N] = find_multiple([3;11;14;7;12;14;8],[(1:10)';(1:10)';3],0)
%				i = [		3		13		21		7		17		8		18]'
%				N = [		3		3		3		7		7		8		8]'
%
%			-> [i,N] = find_multiple([3;11;14;7;12;14;8],[(1:10)';(1:10)';3],1)
%				i = [		3		0		0		7		0		0		8
%    						13		0		0		17		0		0		18
%						   21		0		0		0		0		0		0]
%				N = [3	3	3	7	7	8	8]'
%
%			-> [i,N] = find_multiple([3;11;14;7;12;14;8],[(1:10)';(1:10)';3],2)
%				i = [		3		13		21
%    						7		17		0
% 						   8		18		0]
%				N = [3	7	8]'

% Woodox 2.00 - 31/01/01 - Olivier Vitrac - rev. 08/05/11

% Revision history
% 08/05/11 fix one single entry (with 1 or more occurrence) is found


if nargin<3, mise_en_forme = 0; end

rech =rech(:); base=base(:);

[ind,posc] = find(repmat(rech',length(base),1)==repmat(base,1,length(rech)));
if mise_en_forme && any(ind)
   idx = find([0;diff(posc)]);
   if any(idx)
      idx2 = [idx(1)-1;diff(idx)];
      pas = ones(size(posc));
      pas(idx)=-idx2+1;
      posl = cumsum(pas);
      tab = zeros(max(posl),max(posc));
      i_tab = sub2ind(size(tab),posl,posc);
      tab(i_tab) = ind;
      ind = tab;
   else % added Olivier 08/05/11
      ind = [zeros(size(ind,1),posc(1)-1) ind];
   end
end

if mise_en_forme == 2
   [i,j] = find(ind); %#ok<ASGLU>
   ind = ind(:,unique(j));
end


if nargout>1, if mise_en_forme == 2, val_trouve = rech(unique(posc)); else val_trouve = rech(posc); end, end
if nargout>2, val_non_trouve = setdiff(rech, val_trouve); end