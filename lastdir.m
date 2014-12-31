function [rep,root] = lastdir(chemin)
% LASTDIR extrait le nom du dernier répertoire du chemin (et la racine correspondante)
%		ex. rep = last_dir(chemin)
%		options : [rep,root] = last_dir(chemin)

% Woodox 1.0 - 27/02/01 - Olivier Vitrac

if length(chemin)<2
   warning(['''' chemin 'n''est pas un nom de chemin valide.'])
   rep = '';
end

chemin = remove_slash(chemin);
ind = find(chemin == '\');

if any(ind) & ind(end)<length(chemin)
   rep = chemin(ind(end)+1:end);
   if nargout>1
      root = remove_slash(chemin(1:ind(end)));
   end
else
   rep = chemin;
   if nargout>1
      root = '';
   end
end


function chemin = remove_slash(chemin_sl)
if chemin_sl(end)=='\' & chemin_sl(end-1)~=':' ;
   chemin = chemin_sl(1:end-1);
else
   chemin = chemin_sl;
end