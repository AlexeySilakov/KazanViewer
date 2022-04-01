function outstr = kv_scripteditor(varargin)
% kv_scripteditor Application M-file
%    Intended to use with 'kazan' viewer
% KAZAN dataviewer with plugins By Boris Epel & Alexey Silakov
% MPI of Bioinorganic Chemistry, Muelhaim an der Ruhr, 2003
% Free for non-commercial use. Use the program at your own risk. 
% The authors retains all rights. 
% Contact: epel@mpi-muelheim.mpg.de

if nargin == 0 | nargin == 1 % LAUNCH GUI
    
    handles.parentFiga = figure('Units', 'points', 'Position', [200 100 200 400], 'MenuBar', 'none');
    set(handles.parentFiga, 'Tag', 'parentFiga');
    set(handles.parentFiga, 'Name', 'Script Editor');
    set(handles.parentFiga, 'PaperUnits', 'normalized ');
    set(handles.parentFiga, 'CloseRequestFcn', 'kv_scripteditor(''Close_Callback'', guidata(gcbo))');
    set(handles.parentFiga, 'ResizeFcn', 'kv_scripteditor(''Resize_Callback'', guidata(gcbo))');
    brd = 0.01;
    handles.eList = uicontrol(handles.parentFiga, 'Style', 'edit', 'Units', 'points ', 'Max' , 100);
    set(handles.eList, 'Tag', 'eList');
    set(handles.eList, 'Position', [brd, 2*brd+0.05, 1-2*brd, 1-3*brd-0.05]);
    set(handles.eList, 'Callback', 'kv_scripteditor(''eList_Callback'', guidata(gcbo))');
    set(handles.eList, 'HorizontalAlignment', 'left', 'FontName', 'Monospace', ...
    'FontUnits', 'points', 'FontSize', 12);
    if nargin == 1
        set(handles.eList, 'String', varargin{1});
    else
        set(handles.eList, 'String', '');
    end
    set(handles.eList, 'BackgroundColor', 'w');

    handles.pbOK = uicontrol(handles.parentFiga, 'Style', 'pushbutton', 'Units', 'points');
    set(handles.pbOK, 'Tag', 'pbOK');
    set(handles.pbOK, 'Position', [brd, brd, 0.2, 0.05]);
    set(handles.pbOK, 'Callback', 'kv_scripteditor(''pbOK_Callback'', guidata(gcbo))');
    set(handles.pbOK, 'String', 'OK');
    handles.pbCancel = uicontrol(handles.parentFiga, 'Style', 'pushbutton', 'Units', 'points');
    set(handles.pbCancel, 'Tag', 'pbCancel');
    set(handles.pbCancel, 'Position', [0.2+2*brd, brd, 0.2, 0.05]);
    set(handles.pbCancel, 'Callback', 'kv_scripteditor(''pbCancel_Callback'', guidata(gcbo))');
    set(handles.pbCancel, 'String', 'Cancel'); 

    Resize_Callback(handles);
    handles.String = '';
    set(handles.parentFiga, 'Resize', 'on');
    guidata(handles.parentFiga, handles);
    outstr = '';
    try
        waitfor(handles.parentFiga, 'UserData', 'finish');
        handles = guidata(handles.parentFiga);
        outstr = handles.String;
    end
    try
        delete(handles.parentFiga);    
    end
elseif nargin >1&ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
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

function eList_Callback(handles)

function Resize_Callback(handles)
fp = get(handles.parentFiga, 'Position');
brd = 5;
set(handles.pbOK, 'Position', [brd, brd, 30, 15]);
set(handles.pbCancel, 'Position', [2*brd+30, brd, 30, 15]);
set(handles.eList, 'Position', [brd, 15+2*brd, fp(3)-2*brd, fp(4)-3*brd-15]);

function Close_Callback(handles)
set(handles.parentFiga, 'UserData', 'finish');

function pbOK_Callback(handles)
handles.String = get(handles.eList, 'String');
guidata(handles.parentFiga, handles)
set(handles.parentFiga, 'UserData', 'finish');


function pbCancel_Callback(handles)
set(handles.parentFiga, 'UserData', 'finish');
