% Trace selector GUI for KAZAN viewer.
% Internal routine.
% 
% trace = seltrace(cellarray_of_colors, cellarray_of titles)

% KAZAN dataviewer with plugins By Boris Epel & Alexey Silakov
% MPI of Bioinorganic Chemistry, Muelhaim an der Ruhr, 2003
% Free for non-commercial use. Use the program at your own risk. 
% The authors retains all rights. 
% Contact: epel@mpi-muelheim.mpg.de

%==========================================================================
function trace = seltrace(varargin)
%==========================================================================

if nargin < 3 % launch dialog
  colors = varargin{1};
  nbuttons = size(colors, 1)+1;
  if nargin >= 2, titles = varargin{2}; 
  else, for kk=1:nbuttons,titles{kk, 1}=''; end;
  end;
  if ~nbuttons, trace=0; return; end;
  nbrows = floor((nbuttons-1)/4)+1;
  if nbuttons >= 4, ncols = 4;
  else, ncols = rem(nbuttons, 4);
  end
  bsizex = 70;
  bsizey = 18;
  border = 5;
  fig = dialog('MenuBar', 'none', 'Resize', 'on', ...
    'Position', [200, 200, (bsizex+border)*4+border, (bsizey+border)*nbrows+border],...
    'Name','Select trace');
  for k = 1:nbrows
    for l = 1:ncols
      nb = (k-1)*ncols + l;
      if nb < nbuttons
        h = uicontrol('Position', [border+(l-1)*(bsizex+border), border+(nbrows-k)*(bsizey+border), bsizex, bsizey],...
          'Tag', ['b', num2str(k), num2str(l)], 'Parent', fig,...
          'Style', 'pushbutton', 'String', titles{nb}, 'UserData', nb,...
          'BackgroundColor', colors{nb}, 'HandleVisibility', 'off',...
          'ForegroundColor', [1,1,1]-colors2RGB(colors{nb}), 'FontWeight', 'bold',...
          'Callback', 'seltrace(''button_callback'', gcbo, 0)');
      else
        h = uicontrol('Position', [border+(3)*(bsizex+border), border+(nbrows-k)*(bsizey+border), bsizex, bsizey],...
          'Tag', ['b', num2str(k), num2str(l)], 'Parent', fig,...
          'Style', 'pushbutton', 'UserData', -1,...
          'String', 'Cancel',...
          'Callback', 'seltrace(''cancel_callback'', gcbo, 0)');
      end
    end
  end
  set(fig, 'WindowStyle', 'modal');
  uiwait(fig)
  trace = get(fig, 'UserData');
  delete(fig);  
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
  
  try
    if (nargout)
      [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
    else
      feval(varargin{:}); % FEVAL switchyard
    end
  catch
    disp(lasterr);
  end
end

function handles = cancel_callback(h, p1)
fig = get(h, 'Parent');
uiresume(fig)

function handles = button_callback(h, p1)
fig = get(h, 'Parent');
set(fig, 'UserData', get(h, 'UserData'));
uiresume(fig)

function RGB=colors2RGB(col)
if ischar(col)
    switch col
    case 'g', RGB=[0,1,0]; 
    case 'r', RGB=[1,0,0]; 
    case 'c', RGB=[1,0,1]; 
    case 'm', RGB=[1,1,0]; 
    case 'y', RGB=[0,1,1]; 
    case 'k', RGB=[0,0,0]; 
    otherwise, RGB=[0,0,1];
    end
else
    [s1,s2] = size(col);
    if s1 > s2,  RGB=col';
    else RGB=col;
    end
end
