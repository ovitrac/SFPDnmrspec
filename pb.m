function  p = pb(x,m,s,mode)
%PB passe/coupe bande utilisant une distribution normale cumulée (PNORM)
%		syntaxe : p = PB(x,m,s)
%			x : frequences ou table d'indice
%			m : moyenne(s)
%			s : ecart-type(s)
%		ex. pb(0:100,50,-4) - passe-bas
%		ex. pb(0:100,50,4) - passe-haut
%		ex. pb(0:100,[10 60],4) - passe-bande symétrique
%		ex. pb(0:100,[10 60],4) - coupe-bande symétrique
%		ex. pb(0:100,[10 60],[2 -12]) - passe-bande disymétrique
%		ex. pb(0:150,[10 20 80 100 110],1) fonctionne aussi
%		options :  p = PB(x,m,s,mode)
%		mode = 'p' (défaut) utilise une loi normale cumulée
%		mode = 'd' utilisela fonction de densité de la loi normale

% Woodox 3.00 - 05/04/01 - Olivier Vitrac - rev. 15/07/10

% Revision hostory (for MS 2.1)
% 15/07/10 fix s sizes


if nargin<4, mode = 'p'; end
if nargin<3, s = 1; end
if nargin<2, m = 0; end

nm = numel(m);
ns = numel(s);
if ns<nm
    if ns==1
        s = repmat([s(1);-s(1)],floor(nm/2),1);
        if length(s)<nm, s = [s;s(1)]; end
    else
        s = s(mod((1:nm)-1,ns)+1);
    end
end

flip = sign(s(1)); if flip == 0, flip = 1; end

switch lower(mode)
case 'p'
   p1 = pnorm(x,m(1),abs(s(1)),flip<0);
case 'd'
   p1 = dnorm(x,m(1),abs(s(1)),1,flip<0);
end

if length(s)>1
   p2 = pb(x,m(2:end),s(2:end),mode);
   p = flip*min(flip*[p1(:) p2(:)],[],2);
else
   p = p1(:);
end
