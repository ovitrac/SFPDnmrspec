function [pos,larg,amp] = monotone(X,type,zero)
% MONOTONE recherche les segments monotones dans un vecteur
%		ex. iX = monotone(X)
%		options : [pos,larg,amp] = monotonie(X,type,zero)
%			entrées :
%				X = vecteur colonne
%				type = '+' (défaut), '-' ou '0' ('+-' is also implemented)
%				zero = valeur max de dX assimimée 0
%			sorties :
%				pos = position du segment
%				larg = largeur du segment en indices
%				amplitude de variation du segment

% Thermique 1.0 - 30/04/01 - Olivier Vitrac (source : WoodOx 1.22) - rev. 23/05/13

% Revision History
% 26/07/10 add '+-'
% 28/09/12 change default zero from 1e5*eps to max(10*eps,(max(X)-min(X))/1e6)
% 23/05/13 accept empty zero

if nargin<2, type='+'; end
if nargin<3, zero = []; end
if isempty(zero), zero = max(10*eps,(max(X)-min(X))/1e6); end
if strcmp(type,'+-') || strcmp(type,'-+')
    [pp,lp,ap] = monotone(X,'+',zero);
    [pm,lm,am] = monotone(X,'-',zero);
    pos = [pp;pm];
    if nargout>1, larg = [lp;lm]; end
    if nargout>1, amp = [ap;am]; end
    return
end

X = X(:);
dX=diff(X);
S = [(dX>zero)-(dX<-zero)];
switch upper(type)
case '+', indS = find(S>0);
case '-', indS = find(S<0);
case '0', indS = find(S==0);
end
%dindS = [2;diff(indS)]; % not used 28/09/12
if any(indS)
	il = indS([2;diff(indS)]>1);
	ir = indS([diff(indS);2]>1)+1;
else
   il = []; ir = [];
end

if ~nargout && any(il)
   hold on
   plot(X,'b-'), plot(il,X(il),'ro',ir,X(ir),'ms')
   stem(il,X(ir)-X(il),'gd')
else
   if  nargout>0, pos = il; end
   if  nargout>1, larg = ir-il +1; end
   if  nargout>2, amp = X(ir)-X(il); end
end
