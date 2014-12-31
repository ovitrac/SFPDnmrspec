function resizeax(hs,factor)
%RESIZEAX applies a resize factor to axes to shrink or expand them
% Syntax: resizeax(hs,factor)

%MS 2.1 - 29/09/11 - INRA\Olivier Vitrac - rev. 

% arg check
if nargin<2, error('2 arguments are required'); end
if any(~ishandle(hs)), error('some axes are not valid'); end

% for each axes
for i=1:numel(hs)
    pos = get(hs(i),'position');
    W = pos(3); Wnew = W*factor;
    H = pos(4); Hnew = H*factor;
    set(hs(i),'position',[pos(1)-(Wnew-W)/2 pos(2)-(Hnew-H)/2 Wnew Hnew]);
end


