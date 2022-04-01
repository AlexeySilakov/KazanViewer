function varargout = kv_exportimage(varargin)
% kv_exportimage  Print a figure to a file.
%   kv_exportimage(H) opens GUI for writing the figure with handle H 
%   to the file.
%
%   kv_exportimage (without any parameters) operates on the current figure. 
%
%   kv_exportimage(H, 'nogui') writes the figure with handle H to the file
%   FIG<H>.eps where <H> indicates the (integer) handle value and 'eps' is
%   the extension corresponding to the default output file format
%   (Encapsulated PostScript).
%
%   kv_exportimage(H, 'nogui', PARAM1, VALUE1, PARAM2, VALUE2, ...) 
%       creates output image according to the
%       speceified in 'PAR' settings (no GUI). PAR is a structure. Possible fields and values: 
%        'preview': <0> | 1
%            if [1] creates figure with preview
%        'do_real_job': 0 | <1>
%            if [1] (default) creates output file.
%        'FileType': <'eps'> | 'meta' | 'eps2' | 'jpeg' | 'png' | 'pdf' | 'tiff' |
%           meta_c | 'bitmap_c'
%            format of the output file
%            'eps'      Encapsulated PostScript
%            'eps2'     Encapsulated level 2 Postscript (default)
%            'meta'     Enhanced Meta File (windows only)
%            'jpeg'     JPEG files
%            'png'      Portable Network Graphics
%            'tiff'     Tag Image File Format
%            'PDF'
%            'meta_c'   export image to the clipboard if EMF format
%            'bitmap_c' export image to the clipboard if BMP format
%       'TiffPrev': 0 | <1>
%           specifies includig a preview for EPS files.
%       'FileName':
%           name of the output file. If it is not specified, then file name 
%           'FIG<H>' is used where <H> indicates the (integer) handle value
%       'ColorMode': <'Color'> | 'Black/White' | 'Grey'
%           color scheme for output image: color (default), black&white
%           ('b/w') or grey 
%
%       'PaperPosition': 2 element vector
%               width and height of the output image in units described in
%               PAR.PaperUnits  
%       'PaperUnits': 'points' | 'inches' | 'centimeters' | 'normalized'  
%               units for 'Position'. 'normalized'  Maps the lower-left
%               corner of the page to (0,0) and the upper-right corner to (1,1). 
%               default value is determined by current MatLab settings (see
%               Page Setup menu of the figure) 
%       'Renderer': 'painters', 'zbuffer' or 'opengl'
%               specifies the renderer to use and must be set to one of
%               the values . The default value is determined by the current
%               Matlab settings  (see Page Setup menu of the figure)
%
%       'FontName':  'auto', or string
%               Font Name
%       'FontMode': <'scaled'> | 'points'
%                  specifies mode for changing font size.
%                  if 'points' all font  will be set to equal size.
%       'FontSize': 
%                   font size. Default 
%       'FontAngle': <'normal'> | 'italic'
%       'FontWeight': <'normal'> | 'bold'
%       'FontEncoding': <'latin1'> | 'adobe'
%               specifying the character encoding of the font. The default
%               is the Matlab default (specified in the PRINT command).
%
%       'LineMode': <'scaled'> | 'points'
%                  specifies mode for changing line size.
%                  if 'points' all lines  will be set to equal width
%       'LineWidth': 
%                   line width in units 'LineMode'
%   kv_exportimage(H, 'gui', PARAM1, VALUE1, PARAM2, VALUE2, ...) or
%   kv_exportimage(H, PARAM1, VALUE1, PARAM2, VALUE2, ...)
%       creates GUI with predefined settings, specified in PARAM.../VALUE...
%       in this case parameters 'preview' and 'do_real_job' have no effect.
% Examples:
%       %first, lets create some picutre:
%           H = figure; plot(sin([1:100]/5));
%       %this command creates preview of the image: 
%          kv_exportimage(H, 'nogui', 'preview', 1, 'do_real_job', 0);
%       %following command will creates file with name 'kv_exportimage_demo' in current
%       %MatLab directory:
%          kv_exportimage(H, 'nogui', 'preview', 0, 'do_real_job', 1, ...
%               'FileName', 'kv_exportimage_demo.emf', 'FileType', 'meta');
%       %following command has almost all parameters
%          kv_exportimage(H, 'gui', 'FileType', 'meta', 'FileName', 'kv_exportimage_demo.emf', ...
%               'ColorMode', 'Grey', 'PaperPosition', [400, 600], 'PaperUnits', 'points', ...
%               'Renderer', 'zbuffer', 'FontName', 'Arial', 'FontMode', 'points', ...
%               'FontSize', 16, 'FontAngle', 'italic', 'FontWeight', 'bold', ...
%               'FontEncoding', 'adobe', 'LineMode', 'points', 'LineWidth', 2);       
% Alexey Silakov, MPI-BAC 03.01.2005,
% Based on PTF.M by J.Ploeg

callback = 0;
if nargin>=1
    if ischar(varargin{1})
%        is it callback or just start up
        callback = 1;
    end
end
make_gui = -1;
% manage with input parameters and decide, what to do
if callback % load program
    if ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
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
else % Launch gui / execute program
    par = []; % structure for input parameters (can be empty)
    switch nargin
        case {0, 1}
            make_gui = 1;
        case 2
            if ischar(varargin{2})
                if strcmp(varargin{2}, 'nogui')
                    make_gui = 0;
                elseif strcmp(varargin{2}, 'gui')
                    make_gui = 1;
                else
                    error('In this case second parameter must be either ''gui'' or ''nogui''');
                end
            else
                error('Second parameter must be ''string''');
            end
        otherwise
            
            handles.figa = varargin{1};
            if ~ischar(varargin{2})
                error('Second parameter must be ''string''');
            end
            switch varargin{2}
                case 'gui' % this makes possible to launch GUI with predefined parameters
                    make_gui = 1;
                    start_num = 2;
                case 'nogui'
                    make_gui = 0;
                    start_num = 2;
                otherwise
                    make_gui = 1;
                    start_num = 1;
            end
            if mod(nargin-start_num, 2)&nargin~=start_num
                error('Wrong number of input parameters', '');
            end
            % create structure of the paramers
            for ci = (1:2:(nargin-start_num))+start_num
                par = setfield(par, varargin{ci}, varargin{ci+1});
            end
    end
end

% Do job
if make_gui==1 % LAUNCH GUI
    bgcolor = get(0, 'DefaultUicontrolBackgroundColor');
    handles.parentFiga = figure('Units', 'points', 'Position', [10 10 600 300], 'MenuBar', 'figure', 'Resize', 'on');
    set(handles.parentFiga, 'Tag', 'parentFiga');
    set(handles.parentFiga, 'Color', bgcolor);
    set(handles.parentFiga, 'Name', 'Export Image');
    set(handles.parentFiga, 'PaperUnits', 'points');
    set(handles.parentFiga, 'CloseRequestFcn', 'kv_exportimage(''Close_Callback'', guidata(gcbo))');
    set(handles.parentFiga, 'NumberTitle', 'off');
    guidata(handles.parentFiga, handles);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Predefined settings %%%%%%%%%%%%%%%%%%%%%%%%  
    %     for future: get(findall(0, 'type', 'figure'), 'Name')
    handles.figa = [];
    figastruct = [];
    axesstruct = [];
    if nargin >= 1, 
        handles.figa = varargin{1}; 
    else
        handles.figa = gcf; 
    end
    figastruct = get(handles.figa);
    aaa = findobj(get(handles.figa, 'Children'), 'type', 'axes');
    if ~isempty(aaa)
        axesstruct = get(aaa(1));
    end
    handles.GUI = 1;
    handles.FileType = safeget(par, 'FileType', 'eps'); % type of the output file
    handles.TiffPrev = safeget(par, 'TiffPreview', 0); % =1 means insert 'tiff'-preview to a 'eps'file
    handles.FileName = safeget(par, 'FileName', ''); % Name of the output file
    handles.ColorMode = safeget(par, 'ColorMode', 'Color'); % possible values 'color', 'b/w', 'grey'
    handles.do_real_job = 0; % =0 means to do only preview, =1 means creating output file
    handles.preview = 1;
    % parameters of the operated figure and its axes
    tPos = safeget(figastruct, 'PaperPosition', get(0,'DefaultFigurePaperPosition'));
        handles.PaperPosition = safeget(par, 'PaperPosition', tPos([3, 4]));
    PaperUnits = safeget(figastruct, 'PaperUnits', get(0,'DefaultFigurePaperUnits'));
        handles.PaperUnits = safeget(par, 'PaperUnits', PaperUnits);
    %handles.Render = safeget(figastruct, 'Render', get(0,'DefaultFigureRender'));
    handles.FontName = safeget(par, 'FontName', 'auto');
    handles.FontSize = safeget(par, 'FontSize', 1);
    handles.ColorMode = safeget(par, 'ColorMode', 'Color');
    FontAngle= safeget(axesstruct, 'FontAngle', get(0,'DefaultAxesFontAngle'));
        handles.FontAngle =safeget(par, 'FontAngle', FontAngle);
    FontWeight = safeget(axesstruct, 'FontWeight', get(0,'DefaultAxesFontWeight'));        
        handles.FontWeight = safeget(par, 'FontWeight', FontWeight);
    handles.FontMode = safeget(par, 'FontMode', 'scaled');
    handles.FontEncoding = safeget(par, 'FontEncoding', 'Latin-1');
    handles.FontStyle = safeget(par, 'FontStyle', 'normal');
    handles.LineMode = safeget(par, 'LineMode', 'scaled');
    handles.LineWidth = safeget(par, 'LineWidth', 1);
    handles.firsttime  =1;
    %%%%%%%%%%%%%%%%%%% creating graphic objects %%%%%%%%%%%%%%%%%%%%%
    % Image frame
    % handles.frmImage = uicontrol(handles.parentFiga, 'Style', 'frame', 'Units', 'points', 'Tag', 'frmImage', ...
    %'ForegroundColor',[0 0 0], 'Enable', 'off', 'BackgroundColor', bgcolor);
    handles.axImage = axes('Parent', handles.parentFiga, 'Units', 'points', 'Tag', 'axImage', ...
        'Layer', 'top', 'XTickMode', 'manual', 'XTick', [],'YTick', [], 'Box', 'on', 'XTickLabel', '');
    %     plot(sin(1:100))
    
    % File frame
    callback_str = 'kv_exportimage(''File_Callback'', gcbo, guidata(gcbo))';
    handles.frmFile = uicontrol(handles.parentFiga, 'Style', 'frame', 'Units', 'points', 'Tag', 'frmFile', ...
        'String', 'File', 'HorizontalAlignment', 'left', 'BackgroundColor', bgcolor);
    handles.txtFile = uicontrol(handles.parentFiga, 'Style', 'text', 'Units', 'points', 'Tag', 'txtFile', ...
        'String', 'File', 'HorizontalAlignment', 'left', 'BackgroundColor', bgcolor);
    handles.txtFileFormat = uicontrol(handles.parentFiga, 'Style', 'text', 'Units', 'points', 'Tag', 'txtFileFormat', ...
        'String', 'Format', 'HorizontalAlignment', 'right', 'BackgroundColor', bgcolor);
    handles.txtFileColor = uicontrol(handles.parentFiga, 'Style', 'text', 'Units', 'points', 'Tag', 'txtFileFormat', ...
        'String', 'Mode', 'HorizontalAlignment', 'right', 'BackgroundColor', bgcolor);
    %     handles.txtFileName = uicontrol(handles.parentFiga, 'Style', 'text', 'Units', 'points', 'Tag', 'txtFileName', ...
    %         'String', 'Name', 'HorizontalAlignment', 'right', 'BackgroundColor', bgcolor);
    
    UserData_Format =   {'epsc',     'eps';... % printing device type, extension
                         'epsc2',    'eps';...
                         'meta',    'emf';...
                         'jpeg',    'jpg';...
                         'png',     'png';...   
                         'tiff',    'tiff';...
                         'pdf',     'pdf'; ...
                         'meta_c',  'Clipboard'; ...
                         'bitmap_c','Clipboard'};
    val = 1;
    BrowseEnable = 'on';
    if ~strcmp(handles.FileType, 'eps'),
        val = find(strcmp(UserData_Format(:, 1), handles.FileType));
        if strcmp(UserData_Format{val, 1}, 'meta_c')|strcmp(UserData_Format{val, 1}, 'bitmap_c')
            handles.FileName = 'Clipboard';
            BrowseEnable = 'off';
        end
        if isempty(val), 
            val = 1; 
            disp('kv_exportimage: wrong FileType. FileType = ''eps'' is set.');
        end
    end
    val_1 = find(strcmp({'Black/White', 'Grey', 'Color'}, handles.ColorMode));
    if isempty(val_1), val_1 = 3; end
    handles.pmFileFormat = uicontrol(handles.parentFiga, 'Style', 'popupmenu', 'Units', 'points', 'Tag', 'pmFileFormat', ...
        'String', {'Encapsulated PostScript', 'Encapsulated Level 2 PS', 'Enhanced MetaFile (Win)', ...
            'JPEG', 'Portable Network Graphics', 'tiff', 'pdf','Clipboard Meta', 'Clipboard Bitmap'},...
        'UserData', UserData_Format, 'Value', val, 'HorizontalAlignment', 'center', 'Callback', callback_str);
    handles.pmFileColor = uicontrol(handles.parentFiga, 'Style', 'popupmenu', 'Units', 'points', 'Tag', 'pmFileFormat', ...
        'String', {'Black/White', 'Grey', 'Color'}, 'Value', val_1,'HorizontalAlignment', 'center', 'Callback', callback_str);
    handles.chkFileTiff = uicontrol(handles.parentFiga, 'Style', 'checkbox', 'Units', 'points', 'Tag', 'chkFile', ...
        'String', 'Tiff Preview', 'HorizontalAlignment', 'left', 'BackgroundColor', bgcolor, 'Callback', callback_str, ...
        'Value', handles.TiffPrev);
    handles.pbFileBrowse = uicontrol(handles.parentFiga, 'Style', 'pushbutton', 'Units','points', 'Tag', 'pbFileBrowse', ...
        'String', 'Browse', 'HorizontalAlignment', 'center', 'BackgroundColor', bgcolor, 'Callback', callback_str, ...
        'Enable', BrowseEnable);
    handles.edFileBrowse = uicontrol(handles.parentFiga, 'Style', 'edit', 'Units','points', 'Tag', 'edFileBrowse', ...
        'String', handles.FileName, 'HorizontalAlignment', 'left', 'BackgroundColor', bgcolor, 'Callback', callback_str, ...
        'Enable', BrowseEnable);
    
    % Figure frame
    callback_str = 'kv_exportimage(''Figure_Callback'', gcbo, guidata(gcbo))';
    handles.frmFigure = uicontrol(handles.parentFiga, 'Style', 'frame', 'Units', 'points', 'Tag', 'frmFigure', ...
        'String', 'File', 'HorizontalAlignment', 'left', 'BackgroundColor', bgcolor);    
    handles.txtFigure = uicontrol(handles.parentFiga, 'Style', 'text', 'Units', 'points', 'Tag', 'txtFigure', ...
        'String', 'Figure', 'HorizontalAlignment', 'left', 'BackgroundColor', bgcolor);    
    
    handles.txtFigureWidth = uicontrol(handles.parentFiga, 'Style', 'text', 'Units', 'points', 'Tag', 'txtFigureWidth', ...
        'String', 'Width', 'HorizontalAlignment', 'right', 'BackgroundColor', bgcolor);    
    handles.txtFigureHeight = uicontrol(handles.parentFiga, 'Style', 'text', 'Units', 'points', 'Tag', 'txtFigureHeight', ...
        'String', 'Height', 'HorizontalAlignment', 'right', 'BackgroundColor', bgcolor);    
    handles.txtFigureUnits = uicontrol(handles.parentFiga, 'Style', 'text', 'Units', 'points', 'Tag', 'txtFigureUnits', ...
        'String', 'Units', 'HorizontalAlignment', 'right', 'BackgroundColor', bgcolor);    
    
    allunits = {'inches', 'points', 'centimeters', 'normalized'};
    val = find(strcmp(allunits, handles.PaperUnits));
    if isempty(val), val = 1; end
    handles.edFigureWidth = uicontrol(handles.parentFiga, 'Style', 'edit', 'Units', 'points', 'Tag', 'edFigureWidth', ...
        'String', num2str(handles.PaperPosition(1)), 'HorizontalAlignment', 'left', 'BackgroundColor', bgcolor, 'Callback', callback_str);    
    handles.edFigureHeight = uicontrol(handles.parentFiga, 'Style', 'edit', 'Units', 'points', 'Tag', 'edFigureHeight', ...
        'String', num2str(handles.PaperPosition(2)), 'HorizontalAlignment', 'left', 'BackgroundColor', bgcolor, 'Callback', callback_str);  
    handles.pmFigureUnits = uicontrol(handles.parentFiga, 'Style', 'popupmenu', 'Units', 'points', 'Tag', 'pmFigureUnits', ...
        'String', allunits, 'Value', val, 'HorizontalAlignment', 'left', 'BackgroundColor', bgcolor, 'Callback', callback_str);      
    
    %Fonts
    callback_str = 'kv_exportimage(''Fonts_Callback'', gcbo, guidata(gcbo))';
    handles.frmFonts = uicontrol(handles.parentFiga, 'Style', 'frame', 'Units', 'points', 'Tag', 'frmFonts', ...
        'String', 'Fonts', 'HorizontalAlignment', 'left', 'BackgroundColor', bgcolor);    
    handles.txtFonts = uicontrol(handles.parentFiga, 'Style', 'text', 'Units', 'points', 'Tag', 'txtFonts', ...
        'String', 'Fonts', 'HorizontalAlignment', 'left', 'BackgroundColor', bgcolor);   
    
    handles.txtFontsMode = uicontrol(handles.parentFiga, 'Style', 'text', 'Units', 'points', 'Tag', 'txtFontsMode', ...
        'String', 'Mode', 'HorizontalAlignment', 'right', 'BackgroundColor', bgcolor);     
    handles.txtFontsSize = uicontrol(handles.parentFiga, 'Style', 'text', 'Units', 'points', 'Tag', 'txtFontsSize', ...
        'String', 'Size', 'HorizontalAlignment', 'right', 'BackgroundColor', bgcolor);     
    handles.txtFontsName = uicontrol(handles.parentFiga, 'Style', 'text', 'Units', 'points', 'Tag', 'txtFontsName', ...
        'String', 'Name', 'HorizontalAlignment', 'right', 'BackgroundColor', bgcolor);     
    handles.txtFontsStyle = uicontrol(handles.parentFiga, 'Style', 'text', 'Units', 'points', 'Tag', 'txtFontsStyle', ...
        'String', 'Style', 'HorizontalAlignment', 'right', 'BackgroundColor', bgcolor);     
    handles.txtFontsEncoding = uicontrol(handles.parentFiga, 'Style', 'text', 'Units', 'points', 'Tag', 'txtFontsEncoding', ...
        'String', 'Encoding', 'HorizontalAlignment', 'right', 'BackgroundColor', bgcolor);     
    
    val_1 = find(strcmp({'scaled', 'points'}, handles.FontMode));
    if isempty(val_1), val_1 = 1; end
    val_2 = find(strcmp({'normal', 'bold', 'italic', 'bold/italic'}, handles.FontStyle));
    if isempty(val_2), val_2 = 1; end
    val_3 = find(strcmp({'Latin-1', 'Adobe'}, handles.FontEncoding));
    if isempty(val_3), val_3 = 1; end
    
    handles.pmFontsMode = uicontrol(handles.parentFiga, 'Style', 'popupmenu', 'Units', 'points', 'Tag', 'pmFontsMode', ...
        'String', {'scaled', 'points'}, 'HorizontalAlignment', 'left', 'BackgroundColor', bgcolor, 'Callback', callback_str, ...
        'Value', val_1);     
    handles.edFontsSize = uicontrol(handles.parentFiga, 'Style', 'edit', 'Units', 'points', 'Tag', 'edFontsSize', ...
        'String', num2str(handles.FontSize), 'HorizontalAlignment', 'left', 'BackgroundColor', bgcolor, 'Callback', callback_str);     
    handles.edFontsName = uicontrol(handles.parentFiga, 'Style', 'edit', 'Units', 'points', 'Tag', 'edFontsName', ...
        'String', handles.FontName, 'HorizontalAlignment', 'left', 'BackgroundColor', bgcolor, 'Callback', callback_str);     
    handles.pmFontsStyle = uicontrol(handles.parentFiga, 'Style', 'popupmenu', 'Units', 'points', 'Tag', 'pmFontsStyle', ...
        'String', {'normal', 'bold', 'italic', 'bold/italic'}, 'HorizontalAlignment', 'left', ...
        'BackgroundColor', bgcolor, 'Callback', callback_str, 'Value', val_2);      
    handles.pmFontsEncoding = uicontrol(handles.parentFiga, 'Style', 'popupmenu', 'Units', 'points', 'Tag', 'pmFontsEncoding', ...
        'String', {'Latin-1', 'Adobe'}, 'HorizontalAlignment', 'left', 'BackgroundColor', bgcolor, 'Callback', callback_str, ...
        'Value', val_3);      
    
    %Line
    callback_str = 'kv_exportimage(''Lines_Callback'', gcbo, guidata(gcbo))';
    handles.frmLines = uicontrol(handles.parentFiga, 'Style', 'frame', 'Units', 'points', 'Tag', 'frmLines', ...
        'String', 'Lines', 'HorizontalAlignment', 'left', 'BackgroundColor', bgcolor);    
    handles.txtLines = uicontrol(handles.parentFiga, 'Style', 'text', 'Units', 'points', 'Tag', 'txtLines', ...
        'String', 'Lines', 'HorizontalAlignment', 'left', 'BackgroundColor', bgcolor);   
    
    handles.txtLinesMode = uicontrol(handles.parentFiga, 'Style', 'text', 'Units', 'points', 'Tag', 'txtLinesMode', ...
        'String', 'Mode', 'HorizontalAlignment', 'right', 'BackgroundColor', bgcolor);     
    handles.txtLinesWidth = uicontrol(handles.parentFiga, 'Style', 'text', 'Units', 'points', 'Tag', 'txtLinesWidth', ...
        'String', 'Width', 'HorizontalAlignment', 'right', 'BackgroundColor', bgcolor);    
    
    val_1 = find(strcmp({'scaled', 'points'}, handles.LineMode));
    if isempty(val_1), val_1 = 1; end   
    handles.pmLinesMode = uicontrol(handles.parentFiga, 'Style', 'popupmenu', 'Units', 'points', 'Tag', 'pmLinesMode', ...
        'String', {'scaled', 'points'}, 'Value', val_1, 'HorizontalAlignment', 'left', 'BackgroundColor', bgcolor, ...
        'Callback', callback_str);     
    handles.edLinesWidth = uicontrol(handles.parentFiga, 'Style', 'edit', 'Units', 'points', 'Tag', 'edLinesWidth', ...
        'String', handles.LineWidth, 'HorizontalAlignment', 'left', 'BackgroundColor', bgcolor, 'Callback', callback_str);     
    
    %%%%%%%%%%%%%%%% global checkbox %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    callback_str = 'kv_exportimage(''Store_Callback'', gcbo, guidata(gcbo))';
        handles.chkLeaveChanges = uicontrol(handles.parentFiga, 'Style', 'checkbox', 'Units', 'points', 'Tag', 'chkLeaveChanges', ...
        'String', 'Leave Changes', 'HorizontalAlignment', 'left', 'BackgroundColor', bgcolor, 'Callback', callback_str);  
    % Rendering
    % right now here is nothing
    callback_str = '';%'kv_exportimage(''listFigHands_Callback'', gcbo, guidata(gcbo))';

    AllFigHands = findobj('Type', 'figure');
    idx = find(AllFigHands==handles.figa);
    if isempty(idx), val = 1; else val = idx; end
    str = {};
    for jj = 1:length(AllFigHands)
        
%         str{jj} = [num2str(AllFigHands(jj)), ' ',get(AllFigHands(jj), 'Name')];
        str{jj} = [num2str(get(AllFigHands(jj), 'Number')), ' ',get(AllFigHands(jj), 'Name')];
        
        userdata(jj) = AllFigHands(jj);
    end
    handles.listFigHands = uicontrol(handles.parentFiga, 'Style', 'popupmenu', 'Units', 'points', 'Tag', 'listFigHands', ...
        'String', str, 'Value', val, 'HorizontalAlignment', 'left', 'BackgroundColor', bgcolor, 'Callback', callback_str,...
        'UserData', userdata); 
    handles.txtFigHands = uicontrol(handles.parentFiga, 'Style', 'text', 'Units', 'points', 'Tag', 'txtFigHands', ...
        'String', 'Figure handle:', 'HorizontalAlignment', 'left', 'BackgroundColor', bgcolor);  
    
    % Preview
    callback_str = 'kv_exportimage(''Preview_Callback'', gcbo, guidata(gcbo))';
    handles.frmPrev = uicontrol(handles.parentFiga, 'Style', 'frame', 'Units', 'points', 'Tag', 'frmPrev', ...
        'String', 'Lines', 'HorizontalAlignment', 'left', 'BackgroundColor', bgcolor);    
    handles.txtPrev = uicontrol(handles.parentFiga, 'Style', 'text', 'Units', 'points', 'Tag', 'txtPrev', ...
        'String', 'Preview', 'HorizontalAlignment', 'left', 'BackgroundColor', bgcolor);   
    
    handles.chkPrevAuto = uicontrol(handles.parentFiga, 'Style', 'checkbox', 'Units', 'points', 'Tag', 'chkPrevAuto', ...
        'String', 'Auto preview', 'HorizontalAlignment', 'left', 'BackgroundColor', bgcolor, 'Callback', callback_str);    
    handles.pbPrevPrev = uicontrol(handles.parentFiga, 'Style', 'pushbutton', 'Units', 'points', 'Tag', 'pbPrevPrev', ...
        'String', 'Preview', 'HorizontalAlignment', 'center', 'BackgroundColor', bgcolor, 'Callback', callback_str);   
    handles.pbPrevFullPrev = uicontrol(handles.parentFiga, 'Style', 'pushbutton', 'Units', 'points', 'Tag', 'pbPrevFullPrev', ...
        'String', 'Fullscale preview', 'HorizontalAlignment', 'center', 'BackgroundColor', bgcolor, 'Callback', callback_str);  
    
    % OK/Cancel/Apply
    handles.pbOK = uicontrol(handles.parentFiga, 'Style', 'pushbutton', 'Units', 'points', 'Tag', 'pbOK', ...
        'String', 'OK', 'HorizontalAlignment', 'center', 'BackgroundColor', bgcolor, ...
        'Callback', 'kv_exportimage(''OK_Callback'', gcbo, guidata(gcbo))');   
    handles.pbCancel = uicontrol(handles.parentFiga, 'Style', 'pushbutton', 'Units', 'points', 'Tag', 'pbCancel', ...
        'String', 'Cancel', 'HorizontalAlignment', 'center', 'BackgroundColor', bgcolor, ...
        'Callback', 'kv_exportimage(''Cancel_Callback'', gcbo, guidata(gcbo))');   
    handles.pbApply = uicontrol(handles.parentFiga, 'Style', 'pushbutton', 'Units', 'points', 'Tag', 'pbApply', ...
        'String', 'Apply', 'HorizontalAlignment', 'center', 'BackgroundColor', bgcolor, ...
        'Callback', 'kv_exportimage(''Apply_Callback'', gcbo, guidata(gcbo))');   
    
    set(handles.parentFiga, 'ResizeFcn', 'kv_exportimage(''Figure_ResizeFcn'', guidata(gcbo))');    
    
    guidata(handles.parentFiga, handles);    
    Figure_ResizeFcn(handles)
    
elseif make_gui==0 % no gui
    handles.GUI = 0;
    handles.preview = safeget(par, 'preview', 0);
    handles.do_real_job = safeget(par, 'do_real_job', 1);
    handles.FileType = safeget(par, 'FileType', 'eps'); % type of the output file
    handles.TiffPrev = safeget(par, 'TiffPrev', 0); % =1 means insert 'tiff'-preview to a 'eps'file
    handles.FileName = safeget(par, 'FileName', ['Figure', num2str(handles.figa),'.eps']); % Name of the output file
    handles.ColorMode = safeget(par, 'ColorMode', 'color'); % possible values 'color', 'b/w', 'grey'
    
    % parameters of the operated figure and its axes
    if isfield(par, 'PaperPosition'), handles.PaperPosition = par.PaperPosition; end
    if isfield(par, 'PaperUnits') handles.PaperUnits = par.PaperUnits ; end
    %if isfield(par 'Render') handles.Render = par.Render ; end
    
    handles.FontName = 'auto';
    handles.FontMode = 'scaled';
    handles.FontSize = 1;        
    handles.FontName = 'auto';
    handles.FontMode = 'scaled';
    handles.FontSize = 1;
    if isfield(par, 'FontName'), handles.FontName = par.FontName; end
    if isfield(par, 'FontSize'), handles.FontSize = par.FontSize; end
    if isfield(par, 'FontMode'), handles.FontMode = par.FontMode; end
    if isfield(par, 'FontAngle'), handles.FontAngle = par.FontAngle; end
    if isfield(par, 'FontWeight'), handles.FontWeight = par.FontWeight; end
    if isfield(par, 'FontEncoding'), handles.FontEncoding = par.FontEncoding; end
    
    handles.LineMode = 'scaled';
    handles.LineWidth = 1;
    if isfield(par, 'LineWidth'), handles.LineSize = par.LineWidth; end
    if isfield(par, 'LineMode'), handles.LineMode = par.LineMode; end
    Local_doit(handles);
end

% ------------------------------------------------------------
function Figure_ResizeFcn(handles)
figPos = get(handles.parentFiga, 'Position');
if figPos(3)<432.0000, return; end
if figPos(4)<220.0000, return; end

brd = 5;
lbFrameX = 5; 
textLength = 35;
%%%%%%%%%%%%%%%% Image Frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imFramePos = [brd, brd, (figPos(3)-2*brd)/2, figPos(4)-2*brd];
%set(handles.frmImage, 'Position', imFramePos);
set(handles.axImage, 'Position', [2*brd, 2*brd, (imFramePos(3)-brd), imFramePos(4)-brd]);

%%%%%%%%%%%%%%%% File Frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ff = 75;
fileFramePos =  [brd+sum(imFramePos([1, 3])), figPos(4) - brd-ff, (figPos(3)-4*brd)/4, ff];
pbBrowsePos =   [fileFramePos(1)+brd, fileFramePos(2) + brd, 50, 15];
edFormatPos = [sum(pbBrowsePos([1, 3]))+brd, pbBrowsePos(2),  fileFramePos(3)-3*brd-pbBrowsePos(3), 15];

chkTiffPos =    [fileFramePos(1)+brd + textLength+brd, sum(pbBrowsePos([2, 4])), fileFramePos(3)-3*brd-textLength, 15];
textColorPos =  [fileFramePos(1)+brd, sum(chkTiffPos([2, 4])), textLength, 12];
pmColorPos =   [sum(textColorPos([1, 3]))+brd, sum(chkTiffPos([2, 4])), fileFramePos(3)-textLength - 3*brd, 15];
textFormatPos = [fileFramePos(1)+brd, sum(textColorPos([2, 4]))+brd, textLength, 12];
pmFormatPos =   [sum(textFormatPos([1, 3]))+brd, sum(textColorPos([2, 4]))+brd, fileFramePos(3)-textLength - 3*brd, 15];

set(handles.txtFile, 'Position', [fileFramePos(1)+lbFrameX, sum(fileFramePos([2, 4]))-5, 50, 10]);
set(handles.frmFile, 'Position', fileFramePos);
set(handles.pbFileBrowse, 'Position', pbBrowsePos);
set(handles.chkFileTiff, 'Position', chkTiffPos);
set(handles.txtFileFormat, 'Position', textFormatPos);
set(handles.txtFileColor, 'Position', textColorPos);
set(handles.pmFileFormat, 'Position', pmFormatPos);
set(handles.pmFileColor, 'Position', pmColorPos);
set(handles.edFileBrowse, 'Position', edFormatPos);
%%%%%%%%%%%%%%%% Figure Frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ff2 = 60;
figFramePos =   [brd+sum(imFramePos([1, 3])), fileFramePos(2) - brd-ff2, (figPos(3)-4*brd)/4, ff2];
txtFigureUnits= [figFramePos(1)+brd, figFramePos(2)+brd, textLength, 12];
txtFigureH =    [figFramePos(1)+brd, sum(txtFigureUnits([2, 4]))+brd, textLength, 12];
txtFigureW =    [figFramePos(1)+brd, sum(txtFigureH([2, 4]))+brd, textLength, 12];

edFigureW =     [txtFigureW(1)+textLength+2*brd, txtFigureW(2), figFramePos(3)-textLength - 4*brd, 15];
edFigureH =     [txtFigureH(1)+textLength+2*brd, txtFigureH(2), figFramePos(3)-textLength - 4*brd, 15];
pmFigureUnits = [txtFigureUnits(1)+textLength+2*brd, txtFigureUnits(2), figFramePos(3)-textLength - 4*brd, 15];

set(handles.frmFigure, 'Position', figFramePos);
set(handles.txtFigure, 'Position', [figFramePos(1)+lbFrameX, sum(figFramePos([2, 4]))-5, 50, 10]);
set(handles.txtFigureWidth, 'Position', txtFigureW);
set(handles.txtFigureHeight, 'Position', txtFigureH);
set(handles.txtFigureUnits, 'Position', txtFigureUnits);
set(handles.edFigureWidth, 'Position', edFigureW);
set(handles.edFigureHeight, 'Position', edFigureH);
set(handles.pmFigureUnits, 'Position', pmFigureUnits);

%%%%%%%%%%%%%%% Fonts Frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ff3 = 95;
fntFramePos =   [brd+sum(figFramePos([1, 3])), figPos(4) - brd-ff3, (figPos(3)-4*brd)/4-brd, ff3];

txtFontsEnc =   [fntFramePos(1)+brd, fntFramePos(2)+brd, textLength, 12];
txtFontsSt =    [fntFramePos(1)+brd, sum(txtFontsEnc([2, 4]))+brd, textLength, 12];
txtFontsN =     [fntFramePos(1)+brd, sum(txtFontsSt([2, 4]))+brd, textLength, 12];
txtFontsSi =    [fntFramePos(1)+brd, sum(txtFontsN([2, 4]))+brd, textLength, 12];
txtFontsM =     [fntFramePos(1)+brd, sum(txtFontsSi([2, 4]))+brd, textLength, 12];

pmFontsEnc =    [txtFontsEnc(1)+textLength+2*brd, txtFontsEnc(2), fntFramePos(3)-textLength - 4*brd, 15];
pmFontsSt =     [txtFontsSt(1)+textLength+2*brd, txtFontsSt(2), fntFramePos(3)-textLength - 4*brd, 15];
edFontsN =      [txtFontsN(1)+textLength+2*brd, txtFontsN(2), fntFramePos(3)-textLength - 4*brd, 15];
edFontsSi =     [txtFontsSi(1)+textLength+2*brd, txtFontsSi(2), fntFramePos(3)-textLength - 4*brd, 15];
pmFontsM =      [txtFontsM(1)+textLength+2*brd, txtFontsM(2), fntFramePos(3)-textLength - 4*brd, 15];

set(handles.frmFonts, 'Position', fntFramePos);
set(handles.txtFonts, 'Position', [fntFramePos(1)+lbFrameX, sum(fntFramePos([2, 4]))-5, 50, 10]);
set(handles.txtFontsEncoding, 'Position', txtFontsEnc);
set(handles.txtFontsStyle, 'Position', txtFontsSt);
set(handles.txtFontsName, 'Position', txtFontsN);
set(handles.txtFontsSize, 'Position', txtFontsSi);
set(handles.txtFontsMode, 'Position', txtFontsM);

set(handles.pmFontsEncoding, 'Position', pmFontsEnc);
set(handles.pmFontsStyle, 'Position', pmFontsSt);
set(handles.edFontsName, 'Position', edFontsN);
set(handles.edFontsSize, 'Position', edFontsSi);
set(handles.pmFontsMode, 'Position', pmFontsM);

%%%%%%%%%%%%%% Lines Frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ff4 = 45;
lnsFramePos =   [fntFramePos(1), fntFramePos(2) - ff4-brd, (figPos(3)-4*brd)/4-brd, ff4];

txtLinesWi =    [lnsFramePos(1)+brd, lnsFramePos(2)+brd, textLength, 12];
txtLinesM =     [lnsFramePos(1)+brd, sum(txtLinesWi([2, 4]))+brd, textLength, 12];

edLinesWi =     [txtLinesWi(1)+textLength+2*brd, txtLinesWi(2), lnsFramePos(3)-textLength - 4*brd, 15];
pmLinesM =      [txtLinesM(1)+textLength+2*brd, txtLinesM(2), lnsFramePos(3)-textLength - 4*brd, 15];

set(handles.frmLines, 'Position', lnsFramePos);
set(handles.txtLines, 'Position', [lnsFramePos(1)+lbFrameX, sum(lnsFramePos([2, 4]))-5, 50, 10]);
set(handles.txtLinesWidth, 'Position', txtLinesWi);
set(handles.txtLinesMode, 'Position', txtLinesM);
set(handles.edLinesWidth, 'Position', edLinesWi);
set(handles.pmLinesMode, 'Position', pmLinesM);

set(handles.chkLeaveChanges, 'Position',  [lnsFramePos(1), lnsFramePos(2) - brd-15, (figPos(3)-4*brd)/4-brd, 15]);
set(handles.txtFigHands, 'Position',  [lnsFramePos(1), lnsFramePos(2) - (brd+15)*2, (figPos(3)-4*brd)/4-brd, 15]);
set(handles.listFigHands, 'Position',  [lnsFramePos(1), lnsFramePos(2) - brd*2-15*3, (figPos(3)-4*brd)/4-brd, 15]);

%%%%%%%%%%%% Preveiw Buttons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ff5 = 65;
prvFramePos = [figFramePos(1), figFramePos(2) - brd-ff5, (figPos(3)-4*brd)/4, ff5];
pbPrevFull  = [prvFramePos(1)+brd, prvFramePos(2)+brd, prvFramePos(3)-2*brd, 15];
pbPrevPrev  = [prvFramePos(1)+brd, sum(pbPrevFull([2, 4]))+brd, prvFramePos(3)-2*brd, 15];
chkPrevAuto  = [prvFramePos(1)+brd, sum(pbPrevPrev([2, 4]))+brd, prvFramePos(3)-2*brd, 15];

set(handles.frmPrev, 'Position', prvFramePos);
set(handles.txtPrev, 'Position', [prvFramePos(1)+lbFrameX, sum(prvFramePos([2, 4]))-5, 50, 10]);
set(handles.pbPrevFullPrev, 'Position', pbPrevFull);
set(handles.pbPrevPrev, 'Position', pbPrevPrev);
set(handles.chkPrevAuto, 'Position', chkPrevAuto);

%%%%%%%%%%%%%%%%%%%%%% rest buttons %%%%%%%%%%%%%%%%%%%%%%%%
pbLen = (figPos(3)/2-4*brd)/3;
set(handles.pbApply, 'Position', [figPos(3)-brd-pbLen, brd, pbLen, 20]);
set(handles.pbCancel, 'Position', [figPos(3)-2*brd-2*pbLen, brd, pbLen, 20]);
set(handles.pbOK, 'Position', [figPos(3)-3*brd-3*pbLen, brd, pbLen, 20]);

% ------------------------------------------------------------
function Close_Callback(handles)    
closereq;

% ------------------------------------------------------------
function File_Callback(h, handles)
set(handles.pbFileBrowse, 'Enable', 'on');
set(handles.edFileBrowse, 'Enable', 'on');    
switch h
    case handles.pmFileFormat
        num = get(h, 'Value');
        str = get(h, 'UserData');
        handles.FileType = str{num, 1}; % device
        try
            if strcmp(str{num, 2}, 'Clipboard')
                set(handles.edFileBrowse, 'String', str{num, 2});
                handles.FileName = str{num, 2};
                set(handles.pbFileBrowse, 'Enable', 'off');
                set(handles.edFileBrowse, 'Enable', 'off');                
            else
                [path,name,ext] = fileparts(get(handles.edFileBrowse, 'String'));
                if strcmp(path(end), '\'), 
                    addstr = '';
                else
                    addstr = '\';
                end
                set(handles.edFileBrowse, 'String', [path, addstr, name, '.', str{num, 2}]);
                handles.FileName = [path, addstr, name, '.', str{num, 2}];
            end
        catch
            return;
        end
    case handles.chkFileTiff
        handles.TiffPrev = get(handles.chkFileTiff, 'Value');
    case handles.pbFileBrowse
        [fname,pname] = uiputfile;
        if ~ischar(fname), return; end
        handles.FileName = [pname, fname];
        set(handles.edFileBrowse, 'String', [pname, fname]);
    case handles.edFileBrowse
        handles.FileName = get(h, 'String');
    case handles.pmFileColor
        str = get(h, 'String');
        num = get(h, 'Value');
        handles.ColorMode = str{num};
end
guidata(handles.parentFiga, handles);

% ------------------------------------------------------------
function Figure_Callback(h, handles)
switch h
    case handles.edFigureWidth
        if ~isnumeric(get(h, 'String'))
            handles.PaperPosition(1) = str2num(get(h, 'String'));
        else
            set(h, 'String', '0');
        end
    case handles.edFigureHeight
        if ~isnumeric(get(h, 'String'))
            handles.PaperPosition(2) = str2num(get(h, 'String'));
        else
            set(h, 'String', '0');
        end
    case handles.pmFigureUnits
        str = get(h, 'String');
        num = get(h, 'Value');
        handles.PaperUnits = str{num};
end
guidata(handles.parentFiga, handles);
if get(handles.chkPrevAuto, 'Value'), 
    Local_doit(handles);
end

% ------------------------------------------------------------
function Fonts_Callback(h, handles)
enbl = 1;
switch h
    case handles.pmFontsMode
        str = get(h, 'String');
        num = get(h, 'Value'); 
        handles.FontMode = str{num};
        % for future: follow default settings
        switch str{num}
            case 'scaled'
                set(handles.edFontsSize, 'String', '1');
                handles.FontSize = 1;
            case 'points'
                set(handles.edFontsSize, 'String', '12');
                handles.FontSize = 12;
        end
        enbl = 0;
    case handles.edFontsSize
        if ~isnumeric(get(h, 'String'))
            handles.FontSize = str2num(get(h, 'String'));
        else
            set(h, 'String', '0');
        end
    case handles.edFontsName
        handles.FontName = get(h, 'String');
    case handles.pmFontsStyle
        str = get(h, 'String');
        num = get(h, 'Value');
        switch str{num}
            case 'bold'
                handles.FontWeight = 'bold';
                handles.FontAngle = 'normal';
            case 'italic'
                handles.FontAngle = 'italic';
                handles.FontWeight = 'normal';
            case 'bold/italic'
                handles.FontWeight = 'bold';
                handles.FontAngle = 'italic';
            otherwise
                handles.FontWeight = 'normal';
                handles.FontAngle = 'normal';
        end
    case handles.pmFontsEncoding
end
guidata(handles.parentFiga, handles);
if get(handles.chkPrevAuto, 'Value')&enbl, 
    Local_doit(handles);
end


% ------------------------------------------------------------
function Lines_Callback(h, handles)
switch h
    case handles.edLinesWidth
        if ~isnumeric(get(h, 'String'))
            handles.LineWidth = str2num(get(h, 'String'));
        else
            set(h, 'String', '0');
        end    
    case handles.pmLinesMode
        str = get(h, 'String');
        num = get(h, 'Value'); 
        handles.LineMode = str{num};
end
guidata(handles.parentFiga, handles);
if get(handles.chkPrevAuto, 'Value'), 
    Local_doit(handles);
end

% ------------------------------------------------------------
function Preview_Callback(h, handles)
switch h
    case handles.chkPrevAuto
        if get(handles.chkPrevAuto, 'Value'), 
            Local_doit(handles);
        end
    case handles.pbPrevPrev
        Local_doit(handles);
end
% guidata(handles.parentFiga, handles);

% ------------------------------------------------------------
function OK_Callback(h, handles)
handles.do_real_job = 1; % make real file
guidata(handles.parentFiga, handles);
Local_doit(handles)
Close_Callback(handles);
% ------------------------------------------------------------
function Apply_Callback(h, handles)
handles.do_real_job = 1; % make real file
guidata(handles.parentFiga, handles);
handles = Local_doit(handles);
handles.do_real_job = 1; % return to the preview mode
guidata(handles.parentFiga, handles);
% ------------------------------------------------------------
function Cancel_Callback(h, handles)
Close_Callback(handles);

% ------------------------------------------------------------
function varargout = Local_doit(handles)

num = get(handles.listFigHands, 'Value');
usdata = get(handles.listFigHands, 'UserData');
handles.figa = usdata(num);
if isempty(handles.figa)|~ishandle(handles.figa), 
    disp('Wrong handle. Check the figure');
    return; 
end

allaxes   = findall(handles.figa,'type','axes');
allimages = findall(handles.figa,'type','image');
alllights = findall(handles.figa,'type','light');
alines    = findall(handles.figa,'type','line');
allpatch  = findall(handles.figa,'type','patch');
allrect   = findall(handles.figa,'type','rectangle');
allsurf   = findall(handles.figa,'type','surface');
alltext   = findall(handles.figa,'type','text');
allfont   = [alltext;allaxes];
allcolor  = [alines;alltext;allaxes;alllights];
allmarker = [alines;allpatch;allsurf];
alledge   = [allpatch;allsurf];
allcdata  = [allimages;allpatch;allsurf];
alllines  = [alines; allaxes; allpatch];

% make backup of a properties
figastruct = get(handles.figa);
% old_... {:, 1}=handles, {:, 2}=property name, {:, 3}=value

old_figa = Local_getprop(handles.figa, {'PaperUnits', 'Units', 'Position', 'PaperPosition', 'PaperPositionMode','PaperSize',...
        'PaperType','Render', 'Colormap'});
old_fonts = Local_getprop(allfont, {'FontSize', 'FontName', 'FontAngle', 'FontWeight'});
old.FontSize = old_fonts(1:length(allfont), :); % get only FontSize
old_lines = Local_getprop(alllines, {'LineWidth'});
old_axescoll = Local_getprop(allaxes, {'Color', 'XColor', 'YColor', 'ZColor'});
old_allcol = Local_getprop(allcolor, {'Color'});
old_marker = Local_getprop(allmarker, {'MarkerEdgeColor', 'MarkerFaceColor'}); 
old_edge = Local_getprop(alledge,{'EdgeColor','FaceColor'});
old_cdata = Local_getprop(allcdata,{'CData'});

if handles.firsttime
    handles.old = [old_figa; old_fonts; old_lines; old_axescoll;...
            old_allcol; old_marker; old_edge; old_cdata];
end
% set new properties 
if isfield(handles, 'PaperUnits')
    set(handles.figa, {'PaperUnits', 'Units'}, {handles.PaperUnits, handles.PaperUnits});
end
if isfield(handles, 'PaperPosition')
    set(handles.figa, {'PaperType', 'PaperPositionMode'}, {'a4letter', 'manual'});
    set(handles.figa, 'PaperPosition', [0, 0, handles.PaperPosition]);
   % set(handles.figa, 'Position', [0, 0, handles.PaperPosition]);
    set(handles.figa, 'PaperSize', handles.PaperPosition);
end
%set(handles, 'Render', handles.Render);
if ~strcmp(handles.FontName, 'auto')
    set(allfont, 'FontName', handles.FontName);
end

if strcmp(handles.FontMode, 'scaled')
    if handles.FontSize~=1
        for ci = 1:length(allfont)
            set(allfont(ci), 'FontSize', old.FontSize{ci, 3}*handles.FontSize);
        end
    end
else
    set(allfont, 'FontSize', handles.FontSize);
end
if isfield(handles, 'FontAngle')
    set(allfont, 'FontAngle', handles.FontAngle);
end
if isfield(handles, 'FontWeight')
    set(allfont, 'FontWeight', handles.FontWeight);
end
if strcmp(handles.LineMode, 'scaled')
    if handles.LineWidth~=1
        for ci = 1:length(alllines)
            set(alllines(ci), 'LineWidth', old_lines{ci, 3}*handles.LineWidth);
        end
    end
else
    set(alllines, 'LineWidth', handles.LineWidth);
end

% if it is neccessary, convert all colors to the shades of gray
if strcmpi(handles.ColorMode, 'Grey') 
    oldcmap=get(handles.figa,'Colormap');
    newgrays=0.30*oldcmap(:,1)+0.59*oldcmap(:,2)+0.11*oldcmap(:,3);
    newcmap=[newgrays newgrays newgrays];
    set(handles.figa,'Colormap',newcmap);
    Local_make_it_gray([old_axescoll; old_allcol; old_marker; old_edge; old_cdata]);
elseif strcmpi(handles.ColorMode, 'Black/White')    
    set(alllines, 'Color', 'k')
end
% % this is the one of the posibilities to get image ...
% F = getframe(handles.figa);

FileName = handles.FileName;
if isempty(FileName)
    handles.do_real_job = 0;
end
tfilename = 'kv_temp.png';
tfiletype = 'png';
flags = '';

if strfind(handles.FileType, 'pdf')
    flags = 'write';
end
% Another one

if handles.do_real_job, 
    if strfind(FileName, 'Clipboard'),
        
        print(handles.figa, ['-d', strtok(handles.FileType, '_c')]);
        disp(['Image saved to the clipboard as ', handles.FileType]);
    else
        print(handles.figa, FileName, ['-d',strtok(handles.FileType, '_c'), flags]);
        disp(['File ''', FileName, ''' have been successfully created']);
    end
end
handles.preview = 1;
if handles.preview
    print(handles.figa, tfilename, ['-d', tfiletype]);
    Fpict = imread(tfilename, 'png');
    if isempty(safeget(handles, 'axImage', []))
        figahand = figure;
        handles.axImage = axes('Parent', figahand);
    end
    
    %Fpict = imresize(Fpict,2,'nearest');
    
    image(Fpict, 'Parent', handles.axImage);
    set(handles.axImage, 'XTick', [],'YTick', []);
    delete(tfilename);
    % make equal axis (taken from 'axis')
    a = get(handles.axImage,'Position'); 
    set(handles.axImage,'DataAspectRatio',[1 1 1]);
    % dx = diff(get(handles.axImage,'xlim')); 
    % dy = diff(get(handles.axImage,'ylim'));
    % dz = diff(get(handles.axImage,'zlim'));
    % set(handles.axImage,'PlotBoxAspectRatioMode','auto')
    % pbar = get(handles.axImage,'PlotBoxAspectRatio');
    % set(handles.axImage,'PlotBoxAspectRatio', ...
    %     [a(3) a(4) dz*min(a(3),a(4))/min(dx,dy)]);
end

if ~get(handles.chkLeaveChanges, 'Value')
%     old = [old_figa; old_fonts; old_lines; old_axescoll;...
%         old_allcol; old_marker; old_edge; old_cdata];
    for ci = 1:length(handles.old)
        set(handles.old{ci, 1}, handles.old{ci, 2}, handles.old{ci, 3});
    end
end
handles.firsttime = 0;
guidata(handles.parentFiga, handles);
% figure(handles.parentFiga);
if nargout ==1, varargout{1} = handles; end

% ----------------------------------------------------------
function Local_make_it_gray(input_cells)
% input_cells: {:, 1}= handles, {:, 2}=properties, {:, 3} = old_values
n=length(input_cells(:, 1));
for ci = 1:n
    prop = input_cells{ci, 2};
    if strcmpi(prop,'CData')
        % Map Color Data (Image, Patch or Surface) to shades of gray
        color=input_cells{ci, 3};
        this_is_uint8 = 0;
        if isa(color, 'uint8'),
            % this is about images imported by something like 'imread' 
            color = double(input_cells{ci, 3}) + 1;
            this_is_uint8 = 1;
        end
        if ndims(color)==3 & (isa(color,'double'))
            gray=0.30*color(:,:,1)+0.59*color(:,:,2)+0.11*color(:,:,3);
            color(:,:,1)=gray;
            color(:,:,2)=gray;
            color(:,:,3)=gray;
        end
        if this_is_uint8,
            newvalue = uint8(color - 1);
        else
            newvalue = color;
        end
    else
        % Map color RGB values to shades of gray
        color=input_cells{ci, 3};
        if ~isempty(color)
            if ischar(color)
                switch color(1)
                    case 'y'
                        color=[1 1 0];
                    case 'm'
                        color=[1 0 1];
                    case 'c'
                        color=[0 1 1];
                    case 'r'
                        color=[1 0 0];
                    case 'g'
                        color=[0 1 0];
                    case 'b'
                        color=[0 0 1];
                    case 'w'
                        continue; %color=[1 1 1];
                    case 'k'
                        continue; %color=[0 0 0];
                    case 'auto'
                        continue;
                    otherwise
                        newvalue=color;
                end
            end
            if ~ischar(color)
                color=0.30*color(1)+0.59*color(2)+0.11*color(3);
            end
        end
        if isempty(color) | ischar(color)
            newvalue=color;
        else
            newvalue=[color color color];
        end
    end
    set(input_cells{ci, 1},input_cells{ci, 2},newvalue);
end

% ---------------------------------------------
function outstr = Local_getprop(hand, prop)
outstr ={};
for cj = 1:length(prop)
     for ci = 1:length(hand)
        outstr{end+1, 1} = hand(ci);
        outstr{end, 2} = prop{cj};
        outstr{end, 3} = get(hand(ci), prop{cj});
     end
end

% -----------------------------------------------
function handles = Local_loaddef(handles)
% This function loads default set of the parameters
% parameters of the operated figure and its axes
tPos = safeget(figastruct, 'PaperPosition', get(0,'DefaultFigurePaperPosition'));
tUnits = safeget(figastruct, 'PaperUnits', get(0,'DefaultFigurePaperUnits'));

%handles.Render = safeget(figastruct, 'Render', get(0,'DefaultFigureRender'));
handles.FontName = 'auto';
handles.FontSize = 1;
if isfield(handles, 'FontAngle'), rmfield(handles, 'FontAngle'); end
if isfield(handles, 'FontWeight'), rmfield(handles, 'FontWeight'); end
handles.FontMode = 'scaled';
handles.FontEncoding = 'Latin-1';
handles.LineMode = 'scaled';
handles.LineWidth = 1;
  
set(handles.edFigureWidth, 'String', num2str(tPos(3)));
set(handles.edFigureHeight, 'String', num2str(tPos(4)));

allunits = {'inches', 'points', 'centimeters', 'normalized'};
val = find(strcmp(allunits, tUnits));
set(handles.pmFigureUnits, 'Value', val);
set(handles.pmFontsMode, 'Value', 1);
set(handles.edFontsSize, 'String', num2str(handles.FontSize));
set(handles.edFontsName, 'String', handles.FontName);
set(handles.pmFontsStyle, 'Value', 1);
set(handles.pmFontsEncoding, 'Value', 1);    
set(handles.pmLinesMode, 'Value', 1);
set(handles.edLinesWidth, 'String', '1');
% 

%----------------------------------------------------
function Store_Callback(hObj, handles)