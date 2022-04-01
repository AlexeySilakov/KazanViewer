function varargout = kv_printprofile(varargin)
% KV_PRINTPROFILE utility that stores and applies predefined 
% settings (profiles) for the figures axis, lines and text
% objects.
% Usage:
%   kv_printprofile - launch GUI
%   kv_printprofile(fig_handle, ini_file) - command line
% written for KAZAN dataviewer

% Boris Epel, 2005
% Contact: epel@mpi-muelheim.mpg.de

% Last Modified by GUIDE v2.5 24-Oct-2005 14:01:11

if nargin == 0  % LAUNCH GUI
    
    fig = openfig(mfilename,'reuse');
        set(fig,'Units','pixels');
    
    % Generate a structure of handles to pass to callbacks, and store it. 
    handles = guihandles(fig);
    handles = LoadDefinitions(handles);
    
    handles.FigHandle = 1;
    handles.controls.frame = uicontrol(handles.figure1,...
        'Style', 'frame');
    for ii=1:length(handles.set),
        handles.controls.pagesel(ii) = uicontrol(handles.figure1,...
            'Style', 'togglebutton','String',handles.set{ii}.type,...
            'Callback', 'kv_printprofile(''Page_Callback'',gcbo,[],guidata(gcbo))');
        handles.set{ii}.selected={};
    end
    handles.controls.buttonApply = uicontrol(handles.figure1,...
        'Style', 'pushbutton', 'String', 'Apply',...
        'Callback', 'kv_printprofile(''ApplyOkCancel_Callback'',gcbo,[],guidata(gcbo))');
    handles.controls.buttonOk = uicontrol(handles.figure1,...
        'Style', 'pushbutton', 'String', 'Ok',...
        'Callback', 'kv_printprofile(''ApplyOkCancel_Callback'',gcbo,[],guidata(gcbo))');
    handles.controls.buttonCancel = uicontrol(handles.figure1,...
        'Style', 'pushbutton', 'String', 'Cancel',...
        'Callback', 'kv_printprofile(''ApplyOkCancel_Callback'',gcbo,[],guidata(gcbo))');
    handles.controls.list = uicontrol(handles.figure1,...
        'Style', 'listbox',...
        'Callback', 'kv_printprofile(''List_Callback'',gcbo,[],guidata(gcbo))');
    handles.controls.editnumber = uicontrol(handles.figure1,...
        'Style', 'edit', 'String', '1', ...
        'Callback', 'kv_printprofile(''EditNumber_Callback'',gcbo,[],guidata(gcbo))');
    handles.controls.textnumber = uicontrol(handles.figure1,...
        'Style', 'text', 'String', '1');
    handles.controls.combofld = uicontrol(handles.figure1,...
        'Style', 'popupmenu','String','Nothing');
    handles.controls.buttonAdd = uicontrol(handles.figure1,...
        'Style', 'pushbutton','String','+', ...
        'Callback', 'kv_printprofile(''AddField_Callback'',gcbo,[],guidata(gcbo))');
    handles.controls.buttonDel = uicontrol(handles.figure1,...
        'Style', 'pushbutton','String','-',...
        'Callback', 'kv_printprofile(''DelField_Callback'',gcbo,[],guidata(gcbo))');
    handles.controls.editvalue = uicontrol(handles.figure1,...
        'Style', 'edit','String','','Visible','off',...
        'Callback', 'kv_printprofile(''Value_Callback'',gcbo,[],guidata(gcbo))');
    handles.controls.popupmenuvalue = uicontrol(handles.figure1,...
        'Style', 'popupmenu','String','','Visible','off',...
        'Callback', 'kv_printprofile(''Value_Callback'',gcbo,[],guidata(gcbo))');
    handles.controls.textvtype = uicontrol(handles.figure1,...
        'Style', 'text','String','S');
    
    handles = Page_Callback(handles.controls.pagesel(1), [], handles);
    handles.lastfile.path = '';
    handles.lastfile.name = '';
    handles=figure1_ResizeFcn([], [], handles);
    guidata(fig, handles);
    
    if nargout > 0
        varargout{1} = fig;
    end
    
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
elseif ishandle(varargin{1})  % invoke command line behaviour
    disp(['Running ',varargin{2},' script for Figure ',num2str(varargin{1}),'...']);
    handles = LoadDefinitions([]);
    handles = LoadFile(handles, varargin{2});
    handles.FigHandle = varargin{1};
    ApplyFullSet(handles);
end

%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.


% Code that generates available fields of object 
% figure(1); 

% h = 1; all = fieldnames(get(h));
% h = gca; all = fieldnames(get(h));
% h = plot(1,1); all = fieldnames(get(h));
% all = fieldnames(get(text));
% h = findobj(6,'Type','patch');
% all = fieldnames(get(h(1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% disp('fields_text = {...');
% ii=1;
% while ii < length(all),
%     tline = ''; p = 0;
%     while ii < length(all) & p < 4,
%         if p>0, comma=','; else comma=''; end
%         tline = sprintf('s%s''%s''', tline, comma,all{ii});
%         p = p + 1;
%         ii = ii + 1;
%     end
%     disp([tline,',...']);
% end
% disp('};')    
% 
% disp('fields_type = {...');
% ii=1;
% while ii < length(all),
%     tline = ''; p = 0;
%     while ii < length(all) & p < 4,
%         if p>0, comma=','; else comma=''; end
%         aaa = get(h,all{ii});
%         b = whos('aaa');
%         tline = sprintf('s%s''%s''', tline, comma,b.class);
%         p = p + 1;
%         ii = ii + 1;
%     end
%     disp([tline,',...']);
% end
% disp('};')   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% --------------------------------------------------------------------
function handles = figure1_ResizeFcn(h, eventdata, handles, varargin)

Pos   = get(handles.figure1, 'Position');
Width = Pos(3); Height = Pos(4);
WBorder  = 4;
HBorder  = 4;
HeaderHeight = 20;
HeaderButtonWidth  = 35;
FooterHeight = 30;
FooterButtonWidth  = 60;
FooterButtonHeight  = 25;
EditNumberWidth     = 50;
StandartControlHeight = 20;
VTypeWidth = 10;
set(handles.controls.frame,...
    'Position', [WBorder,HBorder+FooterHeight,Width-2*WBorder,Height-2*HBorder-HeaderHeight+5-FooterHeight]);
for ii=1:length(handles.set),
    set(handles.controls.pagesel(ii),...
        'Position', [WBorder+10+HeaderButtonWidth*(ii-1),Height-HeaderHeight-HBorder,HeaderButtonWidth,HeaderHeight]);
end
set(handles.controls.buttonApply,...
    'Position', [Width-FooterButtonWidth-WBorder,HBorder,FooterButtonWidth,FooterButtonHeight]);
set(handles.controls.buttonOk,...
    'Position', [Width-FooterButtonWidth*3-WBorder*5,HBorder,FooterButtonWidth,FooterButtonHeight]);
set(handles.controls.buttonCancel,...
    'Position', [Width-FooterButtonWidth*2-WBorder*3,HBorder,FooterButtonWidth,FooterButtonHeight]);
set(handles.controls.list,...
    'Position', [2*WBorder,HBorder+2*FooterHeight,Width-4*WBorder,Height-HBorder-4*FooterHeight]);
set(handles.controls.editnumber,...
    'Position', [2*WBorder+1,Height-2*HBorder-HeaderHeight-FooterButtonHeight,EditNumberWidth*2/3,StandartControlHeight]);
set(handles.controls.textnumber,...
    'Position', [2*WBorder+2+EditNumberWidth*2/3,Height-2*HBorder-HeaderHeight-FooterButtonHeight-3,EditNumberWidth/3,StandartControlHeight]);
set(handles.controls.combofld,...
    'Position', [3*WBorder+1+EditNumberWidth,Height-2*HBorder-HeaderHeight-FooterButtonHeight,Width-EditNumberWidth-2*StandartControlHeight-7*WBorder ,StandartControlHeight]);
set(handles.controls.buttonAdd,...
    'Position', [Width-2*StandartControlHeight-3*WBorder-1,Height-2*HBorder-HeaderHeight-FooterButtonHeight,StandartControlHeight,StandartControlHeight]);
set(handles.controls.buttonDel,...
    'Position', [Width-StandartControlHeight-2*WBorder-1,Height-2*HBorder-HeaderHeight-FooterButtonHeight,StandartControlHeight,StandartControlHeight]);
set(handles.controls.editvalue,...
    'Position', [3*WBorder+VTypeWidth,FooterHeight+2*HBorder,Width-5*WBorder-VTypeWidth,StandartControlHeight]);
set(handles.controls.popupmenuvalue,...
    'Position', [3*WBorder+VTypeWidth,FooterHeight+2*HBorder,Width-5*WBorder-VTypeWidth,StandartControlHeight]);
set(handles.controls.textvtype,...
    'Position', [2*WBorder,FooterHeight+2*HBorder-2,VTypeWidth,StandartControlHeight-2]);

% --------------------------------------------------------------------
function handles = LoadDefinitions(handles)
handles.set{5}.fields = {...
        'AlphaDataMapping','AmbientStrength','BackFaceLighting',...
        'CData','CDataMapping',...
        'Clipping',...
        'DiffuseStrength','EdgeAlpha','EdgeColor','EdgeLighting',...
        'EraseMode','FaceAlpha','FaceColor','FaceLighting',...
        'Faces','FaceVertexAlphaData','FaceVertexCData',...
        'LineStyle','LineWidth',...
        'Marker','MarkerEdgeColor','MarkerFaceColor','MarkerSize',...
        'NormalMode',...
        'SpecularColorReflectance','SpecularExponent','SpecularStrength',...
        'VertexNormals',...
        'Vertices','Visible','XData','YData',...
    };
handles.set{5}.fldtypes = {...
        'char','char','char',...
        'double','char',...
        '2',...
        'char','char','char','char',...
        'char','char','char','char',...
        'char','char','char',...
        'char','double',...
        'char','char','char','double',...
        'char',...
        'double','double','double',...
        'double',...
        'double','2','double','double',...
    };
handles.set{5}.type='patch';

handles.set{4}.fields = {...
        'Clipping','Color',...
        'Editing','EraseMode','Extent','FontAngle',...
        'FontName','FontSize','FontUnits','FontWeight',...
        'HorizontalAlignment','Interpreter',...
        'Interruptible','Parent','Position','Rotation',...
        'Selected','SelectionHighlight','String',...
        'Units',...
        'VerticalAlignment',...
    };
handles.set{4}.fldtypes = {...
        '2','double',...
        'char','char','double','3',...
        'char','double','1','4',...
        'char','char',...
        'char','double','double','double',...
        'char','char','char',...
        '1',...
        'char',...
    };
handles.set{4}.type='text';

handles.set{3}.fields = {...
        'Clipping','Color',...
        'EraseMode',...
        'LineStyle','LineWidth','Marker','MarkerEdgeColor',...
        'MarkerFaceColor','MarkerSize','Selected',...
        'SelectionHighlight',...
        'Visible','XData','YData',...
    };
handles.set{3}.fldtypes = {...
        '2','double',...
        'char',...
        'char','double','char','char',...
        'char','double','char',...
        'char',...
        'char','double','double',...
    };
handles.set{3}.type='line';

handles.set{2}.fields = {...
        'ALim','ALimMode','AmbientLightColor',...
        'Box','CameraPosition',...
        'CameraPositionMode','CameraTarget','CameraTargetMode','CameraUpVector',...
        'CameraUpVectorMode','CameraViewAngle','CameraViewAngleMode',...
        'CLim','CLimMode','Clipping','Color',...
        'ColorOrder','DataAspectRatio',...
        'DataAspectRatioMode','DrawMode','FontAngle',...
        'FontName','FontSize','FontUnits','FontWeight',...
        'GridLineStyle',...
        'Layer','LineStyleOrder','LineWidth','MinorGridLineStyle',...
        'NextPlot','PlotBoxAspectRatio','PlotBoxAspectRatioMode',...
        'Position','Projection','Selected','SelectionHighlight',...
        'TickDir','TickDirMode','TickLength',...
        'Units',...
        'View','Visible','XAxisLocation',...
        'XColor','XDir','XGrid','XLabel',...
        'XLim','XLimMode','XMinorGrid','XMinorTick',...
        'XScale','XTick','XTickLabel','XTickLabelMode',...
        'XTickMode','YAxisLocation','YColor','YDir',...
        'YGrid','YLabel','YLim','YLimMode',...
        'YMinorGrid','YMinorTick','YScale','YTick',...
        'YTickLabel','YTickLabelMode','YTickMode','ZColor',...
        'ZDir','ZGrid','ZLabel','ZLim',...
        'ZLimMode','ZMinorGrid','ZMinorTick','ZScale',...
        'ZTick','ZTickLabel','ZTickLabelMode',...
    };
handles.set{2}.fldtypes = {...
        'double','char','double',...
        '2','double',...
        'char','double','char','double',...
        'char','double','char',...
        'double','char','2','double',...
        'double','double',...
        'char','char','3',...
        'char','double','1','4',...
        'char',...
        'char','char','double','char',...
        'char','double','char',...
        'double','char','char','char',...
        '5','char','double',...
        '1',...
        'double','char','char',...
        'double','char','char','double',...
        'double','char','char','char',...
        'char','double','char','char',...
        'char','char','double','char',...
        'char','double','double','char',...
        'char','char','char','double',...
        'char','char','char','double',...
        'char','char','double','double',...
        'char','char','char','char',...
        'double','char','char',...
    };
handles.set{2}.type='axes';

handles.set{1}.fields = {...
        'Alphamap','BackingStore',...
        'Clipping',...
        'Color','Colormap',...
        'CurrentCharacter','CurrentPoint',...
        'Dithermap','DithermapMode','DoubleBuffer','FileName',...
        'FixedColors','IntegerHandle',...
        'InvertHardcopy',...
        'MinColormap','Name','NumberTitle',...
        'PaperOrientation','PaperPosition','PaperPositionMode','PaperSize',...
        'PaperType','PaperUnits','Parent','Pointer',...
        'Position','Renderer',...
        'SelectionHighlight','SelectionType','ShareColors',...
        'Units',...
    };
handles.set{1}.fldtypes = {...
        'double','char',...
        '2',...
        'double','double',...
        'char','double',...
        'double','char','char','char',...
        'double','char',...
        'char',...
        'double','char','char',...
        'char','double','char','double',...
        'char','6','double','char',...
        'double','char',...
        'char','char','char',...
        '1',...
    };
handles.set{1}.type='figure';
handles.options{1}={'inches','centimeters','normalized','points','pixels','characters'};
handles.options{2}={'on','off'};
handles.options{3}={'normal','italic','oblique'};
handles.options{4}={'normal','bold','light','demi'};
handles.options{5}={'in','out'};
handles.options{6}={'normalized','inches','centimeters','points'};

% --------------------------------------------------------------------
function handles = Page_Callback(h, eventdata, handles)
set(handles.controls.pagesel,'Value',0);
set(h,'Value',1);
handles.ActivePage = find(handles.controls.pagesel==h);
set(handles.controls.combofld, 'String',handles.set{handles.ActivePage}.fields, 'Value', 1);
guidata(handles.figure1, handles);
maxval = length(handles.set{handles.ActivePage}.selected);
idx = str2num(get(handles.controls.editnumber, 'String'));
set(handles.controls.editnumber, 'String', num2str(max(min(maxval,idx),1)));
set(handles.controls.textnumber, 'String', num2str(max(1,maxval)));
UpdateList(handles)

% --------------------------------------------------------------------
function varargout = AddField_Callback(h, eventdata, handles)
Str = get(handles.controls.combofld, 'String');
Pos = get(handles.controls.combofld, 'Value');
idx = str2num(get(handles.controls.editnumber, 'String'));

if isempty(handles.set{handles.ActivePage}.selected) |...
        length(handles.set{handles.ActivePage}.selected)<idx,
   handles.set{handles.ActivePage}.selected{idx} = struct(Str{Pos},''); 
elseif ~isfield(handles.set{handles.ActivePage}.selected{idx}, Str{Pos})
    handles.set{handles.ActivePage}.selected{idx} = ...
        setfield(handles.set{handles.ActivePage}.selected{idx}, Str{Pos},'');
end
guidata(handles.figure1, handles);
UpdateList(handles)

% --------------------------------------------------------------------
function varargout = DelField_Callback(h, eventdata, handles)
Pos = get(handles.controls.list,'Value');
idx = str2num(get(handles.controls.editnumber, 'String'));
Str = get(handles.controls.list, 'String');
[fld, str1] = strtok(Str{Pos}, '=');
handles.set{handles.ActivePage}.selected{idx} = ...
    rmfield(handles.set{handles.ActivePage}.selected{idx}, fld);
handles.set{handles.ActivePage}.selected{idx};
guidata(handles.figure1, handles)
UpdateList(handles)

% --------------------------------------------------------------------
function varargout = UpdateList(handles)
Pos = get(handles.controls.list,'Value');
idx = str2num(get(handles.controls.editnumber, 'String'));

if  ~isempty(handles.set{handles.ActivePage}.selected) & ...
        length(handles.set{handles.ActivePage}.selected) >= idx & ...
        isstruct(handles.set{handles.ActivePage}.selected{idx}),
    selfield = fieldnames(handles.set{handles.ActivePage}.selected{idx});
    Str={};
    for ii=1:length(selfield)
        try
            fldidx = getfldidx(handles.set{handles.ActivePage}.fields, selfield{ii});
            Str{end+1}= [selfield{ii},'=', getfield(handles.set{handles.ActivePage}.selected{idx}, selfield{ii})];
        catch
            disp(['Field ''',selfield{ii},''' is wrongly constructed.'])
        end
    end
    set(handles.controls.list, 'String', Str, 'Value', max(1,min(Pos, length(Str))));
else
    set(handles.controls.list, 'String', {});
end

% --------------------------------------------------------------------
function varargout = EditNumber_Callback(h, eventdata, handles)
UpdateList(handles)

% --------------------------------------------------------------------
function varargout = List_Callback(h, eventdata, handles)
idx = str2num(get(handles.controls.editnumber, 'String'));
Str = get(handles.controls.list, 'String');
Pos = get(handles.controls.list,'Value');
[fld, str1] = strtok(Str{Pos}, '=');
val = trim(str1(2:end));
fldidx = getfldidx(handles.set{handles.ActivePage}.fields, fld);
if fldidx > 0,
    type = handles.set{handles.ActivePage}.fldtypes{fldidx};
    switch type
    case 'char',
        set(handles.controls.editvalue,'Visible','on');
        set(handles.controls.popupmenuvalue,'Visible','off');
        set(handles.controls.editvalue,'String', val);
        set(handles.controls.textvtype, 'String','C','TooltipString','String of characters')
    case 'double',
        set(handles.controls.editvalue,'Visible','on');
        set(handles.controls.popupmenuvalue,'Visible','off');
        set(handles.controls.editvalue,'String', val);
        set(handles.controls.textvtype, 'String','N','TooltipString','Numeric value or array')
    otherwise, 
        set(handles.controls.editvalue,'Visible','off');
        Optidx = str2num(type);
        Opt = handles.options{Optidx};
        validx = min(max(1,getfldidx(Opt, val)), length(Opt));
        set(handles.controls.popupmenuvalue,'String', Opt, 'Visible','on', 'Value',validx);
        set(handles.controls.textvtype, 'String','S','TooltipString','Predefined selection');
    end
    set(handles.controls.list, 'UserData', struct('fld',fld,'val',val,'fldidx',fldidx));
end
% getfield(handles.set{handles.ActivePage}.selected{idx},S)

% --------------------------------------------------------------------
function varargout = Value_Callback(h, eventdata, handles)
vvv = get(handles.controls.list, 'UserData');
Str = get(h,'String');
idx = str2num(get(handles.controls.editnumber, 'String'));
switch h
case handles.controls.popupmenuvalue,
    Pos = get(handles.controls.popupmenuvalue,'Value');
    handles.set{handles.ActivePage}.selected{idx} = ...
        setfield(handles.set{handles.ActivePage}.selected{idx},vvv.fld, Str{Pos});
case handles.controls.editvalue,
    handles.set{handles.ActivePage}.selected{idx} = ...
        setfield(handles.set{handles.ActivePage}.selected{idx},vvv.fld, Str);
end
guidata(handles.figure1, handles);
UpdateList(handles);

% --------------------------------------------------------------------
function hh = FindHandleByType(FigHandle, type)
if ~ishandle(FigHandle),
    error('The figure is not specified.');
end
switch type
    case 'figure',
        hh = FigHandle;
    case 'axes',
        hh = findobj(FigHandle,'Type','axes');
    case 'line',
        hhh = findobj(FigHandle,'Type','axes');
        hh  = findobj(hhh, 'Type', 'line');
    case 'text',
        hh  = findobj(FigHandle,'Type','text');
        hhh = findobj(FigHandle,'Type','axes');
        for kk=1:length(hhh)
            hh  = [hh;get(hhh(kk),'Title');get(hhh(kk),'XLabel');get(hhh(kk),'YLabel');get(hhh(kk),'ZLabel')];
        end
    case 'patch'
        hhh = findobj(FigHandle,'Type','axes');
        hh  = findobj(hhh, 'Type', 'patch');
end

% --------------------------------------------------------------------
function varargout = ApplyOkCancel_Callback(h, eventdata, handles)
switch h
    case handles.controls.buttonCancel,
    case {handles.controls.buttonApply, handles.controls.buttonOk}
        ApplyFullSet(handles);
end
if(h==handles.controls.buttonOk | ...
        h==handles.controls.buttonCancel ), closereq; end

% --------------------------------------------------------------------
function ApplyFullSet(handles)
for ii=1:length(handles.set),     % types: figures, axis etc
    if ~isfield(handles.set{ii}, 'selected'), len = 0;
    else len=length(handles.set{ii}.selected);
    end
    hh = FindHandleByType(handles.FigHandle, handles.set{ii}.type);
    for jj=len:-1:1, % sets: 1,2,3 etc
        if jj==1, sethh=hh; 
        else,
            sethh = findobj(hh,'UserData',jj,'Type',handles.set{ii}.type);
            for ll=1:length(sethh)
                hh = hh(hh~=sethh(ll));
            end
        end
        if ~isstruct(handles.set{ii}.selected{jj}), continue; end
        selfield = fieldnames(handles.set{ii}.selected{jj});
        for kk=1:length(selfield)
            fldidx = getfldidx(handles.set{ii}.fields, selfield{kk});
            type = handles.set{ii}.fldtypes{fldidx};
            disp([handles.set{ii}.type, ', set ', num2str(jj), ': ', selfield{kk}]);
            SetValue(sethh, type, selfield{kk}, getfield(handles.set{ii}.selected{jj},selfield{kk}));
        end
    end
end

% --------------------------------------------------------------------
function fldidx = getfldidx(fldlist, fld)
fldidx = -1;
for ii=1:length(fldlist)
    if strcmp(fldlist{ii}, fld),
        fldidx = ii; break;
    end
end

% --------------------------------------------------------------------
function varargout = mFileNew_Callback(h, eventdata, handles, varargin)
for ii=1:length(handles.set),
    handles.set{ii}.selected={};
end
guidata(handles.figure1, handles);
UpdateList(handles);

% --------------------------------------------------------------------
function varargout = mFileLoad_Callback(h, eventdata, handles, varargin)
if isempty(handles.lastfile.path)
    [fname,pname] = uigetfile({'*.ini'}, 'Select ini file name');
else
    curdir = pwd; cd(handles.lastfile.path); 
    [fname,pname] = uigetfile({'*.ini'}, 'Select ini file name');
    cd(curdir);
end

if isequal(fname,0) | isequal(pname,0), return; end
ffull = fullfile(pname, fname);
handles = LoadFile(handles, ffull);
% set(handles.controls.editnumber, 'String', num2str(1));
guidata(handles.figure1, handles);
Page_Callback(handles.controls.pagesel(1), [], handles);

% --------------------------------------------------------------------
function handles = LoadFile(handles, ffull)
ini = inimanaging(ffull);

for ii=1:length(handles.set),
    handles.set{ii}.selected={};
end

inisections = fieldnames(ini);
for ii=1:length(inisections)
    type='';
    for jj=1:length(handles.set),
        if ~isempty(strfind(inisections{ii}, handles.set{jj}.type)),
            type = handles.set{jj}.type; break;
        end
    end
    if isempty(type), continue; end
    if isempty(strfind(inisections{ii},'info'))
        set = str2num(inisections{ii}(length(type)+5:end));
        fld = getfield(ini, inisections{ii});
        handles.set{jj}.selected{set}=fld;
    else
        fld = getfield(ini, inisections{ii});
        handles.set{jj}.info=fld;
    end
end

% --------------------------------------------------------------------
function mFileCopy_Callback(h, eventdata, handles, varargin)
disp('We are cool !');

% --------------------------------------------------------------------
function varargout = mFileSave_Callback(h, eventdata, handles, varargin)
if isempty(handles.lastfile.name), 
    mFileSaveAs_Callback(h, eventdata, handles); return;
end

ini = [];
for ii=1:length(handles.set),
    section = handles.set{ii}.type;
    for jj=1:length(handles.set{ii}.selected)
        item=[];
        sellist = handles.set{ii}.selected{jj};
        if~isempty(sellist)
            sellistfld = fieldnames(sellist);
            for kk=1:length(sellistfld),
                item = setfield(item, sellistfld{kk}, getfield(sellist,  sellistfld{kk}));
            end
            if isempty(item), item.default = 1; end
            ini = setfield(ini, [section,'_set',num2str(jj)], item);
        end
    end
end
for ii=1:length(handles.set)
    if isfield(handles.set{ii}, 'info')
        section = [handles.set{ii}.type,'_info'];
        ini = setfield(ini, section, handles.set{ii}.info);
    end
end
inimanaging(handles.lastfile.name, ini);
disp(['File ',handles.lastfile.name,' saved.']);

% --------------------------------------------------------------------
function varargout = mFileSaveAs_Callback(h, eventdata, handles, varargin)
if isempty(handles.lastfile.path)
    [fname,pname] = uiputfile({'*.ini'}, 'Select ini file name');
else
    curdir = pwd; cd(handles.lastfile.path); 
    [fname,pname] = uiputfile({'*.ini'}, 'Select ini file name');
    cd(curdir);
end
if isequal(fname,0) | isequal(pname,0), return; end
[fpath,name,ext] = fileparts(fname); 
ffull = fullfile(pname, [name, '.ini']);
handles.lastfile.path = pname;
handles.lastfile.name = ffull;
guidata(handles.figure1, handles);

mFileSave_Callback(h, eventdata, handles);
% --------------------------------------------------------------------
function varargout = mHelpHelp_Callback(h, eventdata, handles, varargin)
h = msgbox({'Usage:';'   kv_printprofile - for GUI'; '   kv_printprofile(fig,ini_file) - for command line';'';...
        'Five types of objects are recognised:';...
        'figure, axes, line, text and patch. Axis labels are considered as a text.';'';...
        'Objects are distinguished by type and their ''UserData'' property';...
        '(integer number > 0). For each combination of type and ''UserData'' property';...
        'separate set of parameters can be created.';'';...
        'All objects with ''UserData'' not equal to the number of defined';...
        'set of parameters treated similar to ''UserData''=1 objects';'';...
        ['Axis ''Assign default set ...'' utility assigns following values for label''s ''UserData'' property:',...
            '   title - 2; xlabel - 3; ylabel - 4; zlabel - 5. Correspondingly the parameters of these',...
            ' labels can be edited as a text set with the given above number.'];'';...
        'To add startup menu to any figure, execute the following line:';...
        'uimenu(''Parent'', <fig_number>, ''Label'', ''Profile'', ''Callback'', ''kv_printprofile'')'... 
    },...
    'Help', 'help');
% --------------------------------------------------------------------
function varargout = mHelpAbout_Callback(h, eventdata, handles, varargin)
msgbox({'Very useful utility for some purpose.';'Version unknown.'; 'By Boris Epel 2005.'}, ...
    'About', 'help')
% --------------------------------------------------------------------
function varargout = mFigure_Callback(h, eventdata, handles, varargin)
forms = get(0,'Children');
delete(get(h, 'Children'));

for ii=length(forms):-1:1
    if handles.FigHandle == forms(ii), isactive='on'; else isactive='off'; end    
    uimenu(h, 'Label', ['Figure ', num2str(forms(ii))], 'UserData', forms(ii), ...
        'Checked', isactive, ...
        'Callback', 'kv_printprofile(''SelectFigure_Callback'',gcbo,[],guidata(gcbo))');
end
% --------------------------------------------------------------------
function varargout = SelectFigure_Callback(h, eventdata, handles, varargin)
handles.FigHandle = get(h,'UserData');
disp(['Figure ', num2str(handles.FigHandle), ' activated.'])
guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function varargout = mUtilitiesAxisDefault_Callback(h, eventdata, handles, varargin)
hh  = findobj(handles.FigHandle,'Type','axes');
hhh = get(hh,'Title');
if(~iscell(hhh)) hhh = {hhh}; end
for ii=1:length(hhh)
    set(hhh{ii}, 'UserData',2);
end
hhh = get(hh,'XLabel');
if(~iscell(hhh)) hhh = {hhh}; end
for ii=1:length(hhh)
    set(hhh{ii}, 'UserData',3);
end
hhh = get(hh,'YLabel');
if(~iscell(hhh)) hhh = {hhh}; end
for ii=1:length(hhh)
    set(hhh{ii}, 'UserData',4);
end
hhh = get(hh,'ZLabel');
if(~iscell(hhh)) hhh = {hhh}; end
for ii=1:length(hhh)
    set(hhh{ii}, 'UserData',5);
end

% --------------------------------------------------------------------
function varargout = mUtilitiesLoadFigure_Callback(h, eventdata, handles, varargin)
hh = handles.FigHandle;
ii = 1; jj = 1;
handles.set{ii}.selected{jj} = LoadFields(handles.set{ii}.fields, handles.set{ii}.fldtypes, ...
    hh, handles.set{ii}.selected{jj});
guidata(handles.figure1, handles);
UpdateList(handles);

% --------------------------------------------------------------------
function mUtilitiesLoadAxis_Callback(h, eventdata, handles, varargin)
hhh = handles.FigHandle;
hh  = get(hhh, 'CurrentAxes');
ii = 2; jj = str2num(get(handles.controls.editnumber, 'String'));
handles.set{ii}.selected{jj} = LoadFields(handles.set{ii}.fields, handles.set{ii}.fldtypes, ...
    hh, handles.set{ii}.selected{jj});
guidata(handles.figure1, handles);
UpdateList(handles);

% --------------------------------------------------------------------
function destination = LoadFields(fields, fldtypes, source, destination)
selfield = fieldnames(destination);
for kk=1:length(selfield)
    fldidx = getfldidx(fields, selfield{kk});
    type = fldtypes{fldidx};
    switch type
        case 'double',
            val = get(source,selfield{kk});
            destination = setfield(destination, selfield{kk}, ...
                PrintDouble(val));
        otherwise, 
            destination = setfield(destination, selfield{kk}, ...
                get(source,selfield{kk}));
    end
end

% --------------------------------------------------------------------
function destination = LoadDifferentFields(fields, fldtypes, source1, source2)
destination = [];
for kk=1:length(fields)
    fld = fields{kk};
    val1 = get(source1,fld);
    val2 = get(source2,fld);
    disp(fld)
    if ~ishandle(val1),
        if sum(size(val1)~=size(val2)) | sum(sum(val1~=val2))
            fldidx = getfldidx(fields,fld);
            type = fldtypes{fldidx};
            switch type
                case 'double',
                    destination = setfield(destination, fld, PrintDouble(val2));
                otherwise,
                    destination = setfield(destination, fld, PrintStrings(val2));
            end
        end
    end
end
% --------------------------------------------------------------------
function out = PrintDouble(in)
out = [];
[sz1, sz2]=size(in);
for ii=1:sz1,
    if ii>1, comma=';'; else comma=''; end
    out = [out, comma];
    for jj=1:sz2,
        if jj>1, comma2=','; else comma2=''; end
          out = [out, comma2, num2str(in(ii,jj))];
    end
end

if length(in) > 1, 
    out = ['[',out,']']; 
end
% --------------------------------------------------------------------

function out = PrintStrings(in)
out = [];
[sz1, sz2]=size(in);
if sz1 > 1, apo=''''; else apo=''; end
for ii=1:sz1,
    if ii>1, comma=';'; else comma=''; end
    out = [out, comma, apo,in(ii,:),apo];
end

if sz1 > 1, 
    out = ['[',out,']']; 
end

% --------------------------------------------------------------------
function mUtilGet_Callback(h, eventdata, handles, varargin)
ii = handles.ActivePage;
idx = str2num(get(handles.controls.editnumber, 'String'));
hh = FindHandleByType(handles.FigHandle, handles.set{ii}.type);
if isempty(hh), return; end;

if ~isfield(handles.set{ii}, 'selected'), len = 0;
else len=length(handles.set{ii}.selected);
end

for jj=len:-1:1, % sets: 1,2,3 etc
    if jj==1, sethh=hh; 
    else,
        sethh = findobj(hh,'UserData',jj);
        for ll=1:length(sethh)
            hh = hh(hh~=sethh(ll));
        end
    end
    if jj==idx, break; end
end
sethh = sethh(1);

Str = get(handles.controls.list, 'String');
Pos = get(handles.controls.list,'Value');
[fld, str1] = strtok(Str{Pos}, '=');

val = get(sethh, fld); 
fldidx = getfldidx(handles.set{ii}.fields,fld);
type = handles.set{ii}.fldtypes{fldidx};
switch type
    case 'double',
        handles.set{ii}.selected{idx} = setfield(handles.set{ii}.selected{idx}, fld, PrintDouble(val));
    otherwise,
        handles.set{ii}.selected{idx} = setfield(handles.set{ii}.selected{idx}, fld, PrintStrings(val));
end

guidata(handles.figure1, handles);
UpdateList(handles);
List_Callback(h, eventdata, handles);

% --------------------------------------------------------------------
function mUtilSet_Callback(h, eventdata, handles, varargin)
ii = handles.ActivePage;
idx = str2num(get(handles.controls.editnumber, 'String'));
hh = FindHandleByType(handles.FigHandle, handles.set{ii}.type);
if isempty(hh), return; end;

if ~isfield(handles.set{ii}, 'selected'), len = 0;
else len=length(handles.set{ii}.selected);
end

for jj=len:-1:1, % sets: 1,2,3 etc
    if jj==1, sethh=hh; 
    else,
        sethh = findobj(hh,'UserData',jj);
        for ll=1:length(sethh)
            hh = hh(hh~=sethh(ll));
        end
    end
    if jj==idx, break; end
end

Str = get(handles.controls.list, 'String');
Pos = get(handles.controls.list,'Value');
[fld, str1] = strtok(Str{Pos}, '=');

fldidx = getfldidx(handles.set{ii}.fields, fld);
vtype = handles.set{ii}.fldtypes{fldidx};
disp([handles.set{ii}.type, ', set ', num2str(jj), ': ', fld]);
SetValue(sethh, vtype, fld, getfield(handles.set{ii}.selected{jj},fld));

% --------------------------------------------------------------------
function SetValue(sethh, vtype, fld, value)
switch vtype
    case 'double',
        set(sethh, fld, eval(value));
    otherwise, 
        set(sethh, fld, value);
end

% --------------------------------------------------------------------
function mUtilFindDiff_Callback(h, eventdata, handles, varargin)
hh = figure;
ii=1; jj=1;
handles.set{ii}.selected{jj}=LoadDifferentFields(handles.set{ii}.fields, ...
     handles.set{ii}.fldtypes, hh, handles.FigHandle);
closereq
guidata(handles.figure1, handles);
UpdateList(handles);

% --------------------------------------------------------------------
function mUtilFindDiffAxis_Callback(h, eventdata, handles, varargin)
hhh = figure; hh = axes;
ii=2; jj=str2num(get(handles.controls.editnumber, 'String'));
handles.set{ii}.selected{jj}=LoadDifferentFields(handles.set{ii}.fields, ...
     handles.set{ii}.fldtypes, hh, get(handles.FigHandle, 'CurrentAxes'));
closereq
guidata(handles.figure1, handles);
UpdateList(handles);


% --------------------------------------------------------------------
function mFileExpTiff_Callback(h, eventdata, handles, varargin)
[fname,pname] = uiputfile({'*.tif'}, 'Select file name');
if isequal(fname,0) | isequal(pname,0), return; end
[fpath,name,ext] = fileparts(fname); 
ffull = fullfile(pname, [name, '.tif']);
switch h
    case handles.mFileExpTiff,
        print(['-f',num2str(handles.FigHandle)], '-dtiff', '-r300', ffull)
        disp(['Printing TIFF 300dpi to ',ffull])
    case handles.mFileExpTiff450,
        print(['-f',num2str(handles.FigHandle)], '-dtiff', '-r450', ffull)
        disp(['Printing TIFF 450dpi to ',ffull])
    case handles.mFileExpTiff600,
        print(['-f',num2str(handles.FigHandle)], '-dtiff', '-r600', ffull)
        disp(['Printing TIFF 600dpi to ',ffull])
    case handles.mFileExpTiff900,
        print(['-f',num2str(handles.FigHandle)], '-dtiff', '-r900', ffull)
        disp(['Printing TIFF 900dpi to ',ffull])
end    


% --------------------------------------------------------------------
function mHelpInfo_Callback(hObject, eventdata, handles)
str = {};
for ii=1:length(handles.set)
    str{end+1,1} =  ['   ',upper(handles.set{ii}.type)];
    if isfield(handles.set{ii}, 'info')
        fldnames = fieldnames(handles.set{ii}.info);
        for jj=1:length(fldnames)
            str{end+1,1} = [fldnames{jj}(9:end),' ',getfield(handles.set{ii}.info,fldnames{jj})];
        end
        str{end+1,1} ='';
    end
end
msgbox(str','Info')
