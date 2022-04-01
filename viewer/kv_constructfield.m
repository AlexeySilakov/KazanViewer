function outstr = kv_constructfield(varargin)
% kv_constructfield Application M-file
%    Intended to use with 'kazan' viewer
% KAZAN dataviewer with plugins By Boris Epel & Alexey Silakov
% MPI of Bioinorganic Chemistry, Muelhaim an der Ruhr, 2003
% Free for non-commercial use. Use the program at your own risk. 
% The authors retains all rights. 
% Contact: epel@mpi-muelheim.mpg.de

% Last Modified by GUIDE v2.5 05-Oct-2004 22:40:24
% Alexey Silakov & Boris Epel 5-dec-2003, MPI
% boep 10-dec-2003, MPI

if nargin == 0 | nargin == 1 % LAUNCH GUI
    
    figPos = [500 300 200 100];
    handles.parentFiga = figure('Units', 'pixels', 'Position', figPos, 'MenuBar', 'none');
    set(handles.parentFiga, 'Tag', 'parentFiga');
    set(handles.parentFiga, 'Name', 'Field Create');
    set(handles.parentFiga, 'CloseRequestFcn', 'fittingbox_func(''Close_Callback'', gcbo)');
    
    sp = 5;
    butH = (figPos(4)-4*sp)/3;
    
    handles.eField = uicontrol(handles.parentFiga, 'Style', 'edit', 'Units', 'pixels');
    set(handles.eField, 'Tag', 'eField');
    set(handles.eField, 'Position', [5, 2*butH+3*sp, 190, butH]);
    set(handles.eField, 'Callback', 'kv_constructfield(''eField_Callback'', guidata(gcbo))')
    if nargin == 1 & ischar(varargin{1})
        set(handles.eField, 'String', varargin{1});
    else
        set(handles.eField, 'String', 'Sys.S = 0.5');
    end
    set(handles.eField, 'BackgroundColor', 'w');
    
    handles.popStr = uicontrol(handles.parentFiga, 'Style', 'popupmenu', 'Units', 'pixels');
    set(handles.popStr, 'Tag', 'popStr');
    set(handles.popStr, 'Position', [5, butH+2*sp, 70, butH]);
    set(handles.popStr, 'Callback', 'kv_constructfield(''popStr_Callback'', guidata(gcbo))')
    handles.Str = {'Sys', 'Exp', 'Opt', 'Shift'};
    set(handles.popStr, 'String', handles.Str);
    handles.popPar = uicontrol(handles.parentFiga, 'Style', 'popupmenu', 'Units', 'pixels');
    set(handles.popPar, 'Tag', 'popPar');
    set(handles.popPar, 'Position', [80, butH+2*sp, 115, butH]);
    set(handles.popPar, 'Callback', 'kv_constructfield(''popPar_Callback'', guidata(gcbo))')
    
    handles.pbOK = uicontrol(handles.parentFiga, 'Style', 'pushbutton', 'Units', 'pixels');
    set(handles.pbOK, 'Tag', 'pbOK');
    set(handles.pbOK, 'Position', [5, sp, 60, butH]);
    set(handles.pbOK, 'Callback', 'kv_constructfield(''pbOK_Callback'', guidata(gcbo))')
    set(handles.pbOK, 'String', 'OK');
    handles.pbCancel = uicontrol(handles.parentFiga, 'Style', 'pushbutton', 'Units', 'pixels');
    set(handles.pbCancel, 'Tag', 'pbCancel');
    set(handles.pbCancel, 'Position', [70, sp, 60, butH]);
    set(handles.pbCancel, 'Callback', 'kv_constructfield(''pbCancel_Callback'', guidata(gcbo))')
    set(handles.pbCancel, 'String', 'Cancel'); 
    
    handles.Sysstr = {'A', 'AStrain', 'Afull', 'Apa', 'Bqk', 'D', 'DStrain', 'Dpa', ...
            'HStrain', 'I', 'Nucs', 'n', 'Q', 'Qpa', 'S', 'ee', 'eepa', 'g', ...
            'gStrain', 'gn', 'gpa', 'lw', 'lwEndor'};
    
    handles.Expstr = {'Detection', 'ExciteWidth', 'Field', 'Harmonic', 'Orientation', 'Range',...
            'Temperature', 'gField', 'mwFreq', 'nPoints', 'asym'};
    
    handles.Optstr = {'Enhancement', 'Freq2Field', 'Intensity', 'LineShape', 'Nuclei',...
            'Output', 'PhiRange', 'Sim', 'Symmetry', 'ShowCor', 'ThetaRange', 'Transitions',...
            'Treshold', 'Verbosity', 'nKnots', 'nSpline', 'nTransitions', 'unit_cnv'};
    
    set(handles.popPar, 'String', handles.Sysstr);
    handles.Shiftstr = {'x', 'dx', 'y', 'dy', 's', 'Scale', 'LineStyle', 'LineWidth', 'Marker', 'contour'};
    handles.String = '';
    set(handles.parentFiga, 'Resize', 'off');
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

function eField_Callback(handles)

function popStr_Callback(handles)
num = get(handles.popStr, 'Value');
titlestrs = get(handles.popStr, 'String');
Str = getfield(handles, [handles.Str{num}, 'str']);
set(handles.popPar, 'String', Str);

Fstr = get(handles.eField, 'String');
[names, val] = strtok(Fstr, '=');
[strtype, strname] = strtok(names,'.');
resstr = [titlestrs{num}, strname, val];
set(handles.eField, 'String', resstr);
guidata(handles.parentFiga, handles);
% popPar_Callback(handles)

function popPar_Callback(handles)
num = get(handles.popPar, 'Value');
titlestrs = get(handles.popPar, 'String');
Fstr = get(handles.eField, 'String');
[names, val] = strtok(Fstr, '=');
[strtype, strname] = strtok(names,'.');

resstr = [strtype, '.', titlestrs{num},' ', val];
set(handles.eField, 'String', resstr);
guidata(handles.parentFiga, handles);

function pbOK_Callback(handles)
handles.String = get(handles.eField, 'String');
set(handles.parentFiga, 'UserData', 'finish');
guidata(handles.parentFiga, handles)

function pbCancel_Callback(handles)
set(handles.parentFiga, 'UserData', 'finish');
