% Trace selector GUI for KAZAN viewer.
% Internal routine.
%
% trace = seltrace(cell_array_of_colors, cell_array_of titles)

% KAZAN dataviewer with plugins By Boris Epel & Alexey Silakov
% MPI of Bioinorganic Chemistry, Muelhaim an der Ruhr, 2006
% Free for non-commercial use. Use the program at your own risk.
% The authors retains all rights.
% Contact: epel@mpi-muelheim.mpg.de

%==========================================================================
function trace = kv_seltrace(varargin)
%==========================================================================

if nargin < 2
  warning('kv_seltrace: function is called with less than 2 arguments.');
  return;
elseif nargin < 3 % launch dialog
  colors = varargin{1};
  nbuttons = size(colors, 1)+1;
  if nargin >= 2, titles = varargin{2};
  else for kk=1:nbuttons,titles{kk, 1}=''; end;
  end;
  if ~nbuttons, trace=0; return; end;
  nbrows = floor((nbuttons-1)/4)+1;
  if nbuttons >= 4, ncols = 4;
  else ncols = rem(nbuttons, 4);
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
        uicontrol('Position', [border+(l-1)*(bsizex+border), border+(nbrows-k)*(bsizey+border), bsizex, bsizey],...
          'Tag', ['b', num2str(k), num2str(l)], 'Parent', fig,...
          'Style', 'pushbutton', 'String', titles{nb}, 'UserData', nb,...
          'ForegroundColor', colors{nb}, 'HandleVisibility', 'off',...
          'FontWeight', 'bold',...
          'Callback', 'kv_seltrace(''button_callback'', gcbo, 0)');
      else
        uicontrol('Position', [border+(3)*(bsizex+border), border+(nbrows-k)*(bsizey+border), bsizex, bsizey],...
          'Tag', ['b', num2str(k), num2str(l)], 'Parent', fig,...
          'Style', 'pushbutton', 'UserData', -1,...
          'String', 'Cancel',...
          'Callback', 'kv_seltrace(''cancel_callback'', gcbo, 0)');
      end
    end
  end
  set(fig, 'WindowStyle', 'modal');
  uiwait(fig)
  trace = get(fig, 'UserData');
  delete(fig);
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

  try
      feval(varargin{:}); % FEVAL switchyard
  catch
    disp(lasterr);
  end
end

function cancel_callback(h, p1)
fig = get(h, 'Parent');
uiresume(fig)

function button_callback(h, p1)
fig = get(h, 'Parent');
set(fig, 'UserData', get(h, 'UserData'));
uiresume(fig)

