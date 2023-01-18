function fancycamera(varargin)
%FANCYCAMERA Interactively orbit/zoom/pan with the mouse buttons.

hfig=gcbf;
if isempty(hfig)
  hfig=gcf;
end
if nargin==0
  down=get(hfig, 'windowbuttondownfcn');
  if ~strcmp(down,'fancycamera(''down''  )')
    cameramenu
  end
  set(hfig, 'windowbuttondownfcn',   'fancycamera(''down''  )')
  set(hfig, 'windowbuttonmotionfcn', 'fancycamera(''motion'')')
  set(hfig, 'windowbuttonupfcn',     'fancycamera(''up''    )')
else
  h = findobj(get(hfig,'children'), 'type', 'uimenu', 'tag', 'cm598');
  if isempty(h)
    fancycamera;
    return;
  end
  Udata = get(h(1), 'userdata');
  switch get(hfig,'selectiontype')
   case 'normal'
    Udata.mode='orbit';
   case 'extend'
    Udata.mode='zoom';
   case 'alt'
    Udata.mode='pan';
  end
  set(h(1), 'userdata', Udata);
  cameramenu(varargin{:});
end
