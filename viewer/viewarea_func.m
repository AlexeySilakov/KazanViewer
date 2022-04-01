function viewarea_func(varargin)
% fittingbox_func Application M-file ONLY for fittingbox.m
% Construct GUI for selecting functions
%
% input argument: handles of fittingbox uicontrols
%
% Last Modified 23.02.2004
% boep 23-feb-2004

if nargin == 0
elseif nargin == 1, % launch dialog
    hhh = varargin{1};
    hand = guidata(hhh.MainFigure);

    dlgViewDialog = dialog('Resize', 'on');
    handles = guidata(dlgViewDialog);
    handles.oldraxis = axis(hand.aReal);
    handles.oldiaxis = axis(hand.aImag);
    
    set(dlgViewDialog, 'Tag', 'ViewDialog');
    set(dlgViewDialog, 'Name', 'Select View Area');
    set(dlgViewDialog, 'PaperUnits', 'normalized');
    set(dlgViewDialog, 'Position', [103 200 200 150]);
    set(dlgViewDialog, 'UserData', hhh);
    set(dlgViewDialog, 'CloseRequestFcn', 'viewarea_func(''Close_Callback'', gcbo)');
    

    psize = 0.45;
    brd = 0.04;
    % panels
    tags = {'fImag'; 'fReal'};
    pos(1, :) = [brd, 1-1.5*psize-2*brd, 1-2*brd, psize];
    pos(2, :) = pos(1, :) + [0, psize+brd, 0, -psize/2];
    for k = 1:size(tags, 1)
        h = uicontrol('Parent', dlgViewDialog, 'Tag', tags{k});
        set(h, 'Style', 'frame');
        set(h, 'Units', 'normalized');
        set(h, 'Position', pos(k, :));
    end
    
    % text
    tags = {'tY'; 'tX'; 'tIm'; 'tRe'};
    texts = {'y'; 'x'; 'im'; 're/mag'};
    pos(1, :) = [brd*1.5, 1-psize/2-3*brd-0.02, 0.2, 0.1];
    pos(2, :) = pos(1, :) + [0, psize/2+brd, 0, 0];
    pos(3, :) = [brd*1.1, 1-1.5*psize-brd+0.02, 0.2, 0.1];
    pos(4, :) = pos(3, :) + [0, (psize-2*brd)/2, 0, 0];
    for k = 1:size(tags, 1)
        h = uicontrol('Parent', dlgViewDialog, 'Tag', tags{k});
        set(h, 'Style', 'text');
        set(h, 'Units', 'normalized');
        set(h, 'Position', pos(k, :));
        set(h, 'String', texts{k});
    end
    
    % edit
    tags = {'eIYmin'; 'eRYmin'; 'eRXmin'; ...
            'eIYmax'; 'eRYmax'; 'eRXmax'};
    values = {handles.oldiaxis(3); handles.oldraxis(3); handles.oldraxis(1); 
        handles.oldiaxis(4); handles.oldraxis(4); handles.oldraxis(2)};
    ss = .2 + brd*2;
    sz = (1-brd*5-.2)/2;
    pos(1, :) = [ss, 1-1.5*psize-brd+0.01, sz, 0.12];
    pos(2, :) = pos(1, :) + [0, (psize-2*brd)/2, 0, 0];
    pos(3, :) = pos(1, :) + [0, psize+brd, 0, 0];
    pos(4, :) = pos(1, :) + [sz + brd, 0, 0, 0];
    pos(5, :) = pos(4, :) + [0, (psize-2*brd)/2, 0, 0];
    pos(6, :) = pos(4, :) + [0, psize+brd, 0, 0];
    for k = 1:size(tags, 1)
        h = uicontrol('Parent', dlgViewDialog, 'Tag', tags{k});
        set(h, 'Style', 'edit');
        set(h, 'Units', 'normalized');
        set(h, 'Position', pos(k, :));
        set(h, 'String', num2str(values{k}));
    end
    
    % buttons
    tags = {'pbOk'; 'pbCancel'; 'pbApply'};
    texts = {'Ok'; 'Cancel'; 'Apply'};
    pos(1, :) = [brd,  brd, (1-6*brd)/3, 0.16];
    pos(2, :) = pos(1, :) + [(1-6*brd)/3+2*brd, 0, 0, 0];
    pos(3, :) = pos(2, :) + [(1-6*brd)/3+2*brd, 0, 0, 0];
    for k = 1:size(tags, 1)
        h = uicontrol('Parent', dlgViewDialog, 'Tag', tags{k});
        set(h, 'Style', 'pushbutton');
        set(h, 'Units', 'normalized');
        set(h, 'Position', pos(k, :));
        set(h, 'String', texts{k});
        set(h, 'Callback', ['viewarea_func(''', tags{k}, '_Callback'', gcbo)']);
    end
    
    guidata(dlgViewDialog, handles);
    
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

%---------------------------------------------------
function pbOk_Callback(h)
fighand = get(h, 'Parent');
pbApply_Callback(h);
delete(fighand);

%---------------------------------------------------
function pbCancel_Callback(h)
% Restore previous parameters
fighand = get(h, 'Parent');
handles = guidata(fighand);
hhh = get(fighand, 'UserData');
hand = guidata(hhh.MainFigure);
axis(hand.aReal, handles.oldraxis)
axis(hand.aImag, handles.oldiaxis)
delete(fighand);

%---------------------------------------------------
function pbApply_Callback(h)
fighand = get(h, 'Parent');
handles = guidata(fighand);
hhh = get(fighand, 'UserData');
hand = guidata(hhh.MainFigure);
raxis    = [str2num(get(findobj(fighand, 'Tag', 'eRXmin'),'String')),
    str2num(get(findobj(fighand, 'Tag', 'eRXmax'),'String')),
    str2num(get(findobj(fighand, 'Tag', 'eRYmin'),'String')),
    str2num(get(findobj(fighand, 'Tag', 'eRYmax'),'String'))];
iaxis    = [str2num(get(findobj(fighand, 'Tag', 'eRXmin'),'String')),
    str2num(get(findobj(fighand, 'Tag', 'eRXmax'),'String')),
    str2num(get(findobj(fighand, 'Tag', 'eIYmin'),'String')),
    str2num(get(findobj(fighand, 'Tag', 'eIYmax'),'String'))];
axis(hand.aReal, raxis);
axis(hand.aImag, iaxis);
uiresume

%---------------------------------------------------
function Close_Callback(h)
try
    hand = guidata(get(h, 'UserData'));
    uiresume(hand.figure1);
end
delete(h);