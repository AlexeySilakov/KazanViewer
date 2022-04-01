function varargout = kazan(varargin)
% KAZAN dataviewer with plugins By Boris Epel & Alexey Silakov
% MPI of Bioinorganic Chemistry, Muelhaim an der Ruhr, 2003-2013
% University of Chicago, 2006
% Free for non-commercial use. Use the program at your own risk.
% The authors retain all rights.
% Contact: bepel@uchicago.edu, silakov@mpi-muelheim.mpg.de

% alsi 28-apr-13

% parameters of main form handles
% 'ax.x' column of data x values
% 'ax.xlabel' data x axis unit
% 'y' columns of data y values
% 'rx' column of result x values
% 'ry' columns of result y values
% 'rxlabel' result x axis unit
% 'script' final script that toolbox executes
% 'lscript' loading script

if nargin == 0  % LAUNCH GUI
    oldfig = getappdata(0, 'Kazan_Viewer_open');
    oldfig = oldfig(ishandle(oldfig));
    tver = version;
    MatlabVer = str2num(tver(1:3));
    %     opc = 1; % only one is allowed
    opc = 2; % multiple but no plugins
    %     opc = 3; %

    fpath = fileparts(which('kazan'));
    inifilename = [fpath, filesep, 'kazan.ini'];
    ini = inimanaging(inifilename);
    try opc = str2num(ini.KazanViewer.MultyKazan); catch opc = 1; end
    if opc == 1 && ~isempty(oldfig), set(oldfig,'Visible','on'); return; end
    if isempty(oldfig), opc = 3; end % In order to get plugins

    fig = openfig(mfilename,'new');
    setappdata(0, 'Kazan_Viewer_open', [oldfig, fig]);

    check_figure_size(fig);
    
    % Generate a structure of handles to pass to callbacks, and store it.
    handles = guihandles(fig);

    filename = get(fig, 'FileName');
    [fpath,name,ext] = fileparts(filename);
    handles.inifile = [fpath, filesep, 'kazan.ini'];
    handles.pluginsdir = fpath;
    ini = inimanaging(handles.inifile);

    try handles.dir_path = ini.KazanViewer.lastDir; catch handles.dir_path = pwd; end
    try handles.autoZoomOff = ini.KazanViewer.autoZoomOff; catch handles.autoZoomOff = 'off'; end
    try handles.FixedZoomOpt = ini.KazanViewer.FixedZoomOpt; catch handles.FixedZoomOpt = 'off'; end

    %**********************************************************************
    % Main Menu
    %**********************************************************************
    % menu File
    handles.mFile = uimenu(handles.MainFigure, 'Label','File', 'Tag', 'mFile');
    handles.mFileNewWindow = uimenu(handles.mFile, 'Label','New Window', 'Tag', 'mFileNewWindow', ...
        'Callback', 'kazan(''mFileSave'',gcbo,[],guidata(gcbo))');
    handles.mFileOpen = uimenu(handles.mFile, 'Label','Open Directory', 'Tag', 'mFileOpen', ...
        'Callback', 'kazan(''mFileSave'',gcbo,[],guidata(gcbo))');
    handles.mFileSave = uimenu(handles.mFile, 'Label','Save Source', 'Tag', 'mFileSave', ...
        'Callback', 'kazan(''mFileSave'',gcbo,[],guidata(gcbo))');
    handles.mFileSaveProc = uimenu(handles.mFile, 'Label','Save Source with Processing', 'Tag', 'mFileSaveProc', ...
        'Callback', 'kazan(''mFileSave'',gcbo,[],guidata(gcbo))');
    handles.mFileConvert = uimenu(handles.mFile, 'Label','Convert directory to ASCII', 'Tag', 'mFileConvert', ...
        'Callback', 'kazan(''mFileSave'',gcbo,[],guidata(gcbo))');
    handles.mFileSaveWorksp = uimenu(handles.mFile, 'Label','Export to Workspace', 'Tag', 'mFileSaveWorksp', ...
        'Separator', 'on', 'Callback', 'kazan(''mFileWorkspace'',gcbo,[],guidata(gcbo))');
    handles.mFileLoadWorksp = uimenu(handles.mFile, 'Label','Load from Workspace', 'Tag', 'mFileLoadWorksp', ...
        'Callback', 'kazan(''mFileWorkspace'',gcbo,[],guidata(gcbo))');
    handles.mFavorites = uimenu(handles.mFile, 'Label','Favorite Directories', 'Tag', 'mFavorites', ...
        'Separator', 'on');
    handles.mFavorAdd = uimenu(handles.mFavorites, 'Label','Add', 'Tag', 'mFavorAdd', ...
        'Callback', 'kazan(''mFavorites'',gcbo,[],guidata(gcbo))');
    handles.mClearAll = uimenu(handles.mFavorites, 'Label','Clear All', 'Tag', 'mClearAll', ...
        'Callback', 'kazan(''mFavorites'',gcbo,[],guidata(gcbo))');
    handles.mFile2Source = uimenu(handles.mFile, 'Label','Move to source (<<)', 'Tag', 'mFile2Source', ...
        'Separator', 'on', 'Callback', 'kazan(''mFile2Source'',gcbo,[],guidata(gcbo))');
    handles.mFileAddHistory = uimenu(handles.mFile, 'Label','Add to Selected', 'Tag', 'mFileAddHistory', ...
        'Separator', 'on', 'Callback', 'kazan(''mFileHistory'',gcbo,[],guidata(gcbo))');
    handles.mFileRemoveHistory = uimenu(handles.mFile, 'Label','Remove from Selected', 'Tag', 'mFileRemoveHistory', ...
        'Callback', 'kazan(''mFileHistory'',gcbo,[],guidata(gcbo))');
    handles.mFileRemoveAllHistory = uimenu(handles.mFile, 'Label','Clear Selected files', 'Tag', 'mFileRemoveAllHistory', ...
        'Callback', 'kazan(''mFileHistory'',gcbo,[],guidata(gcbo))');
    handles.mFilePreviousHistory = uimenu(handles.mFile, 'Label','Previous file', 'Tag', 'mFilePreviousHistory', ...
        'Accelerator', 'P', ...
        'Callback', 'kazan(''mFileHistory'',gcbo,[],guidata(gcbo))');
    handles.mFileNextHistory = uimenu(handles.mFile, 'Label','Next file', 'Tag', 'mFileNextHistory', ...
         'Accelerator', 'N', ...
         'Callback', 'kazan(''mFileHistory'',gcbo,[],guidata(gcbo))');
    handles = CreateFavorites(handles);
    % menu View
    handles.mView = uimenu(handles.MainFigure, 'Label','View', 'Tag', 'mView');
    handles.mViewSource = uimenu(handles.mView, 'Label','Source', 'Tag', 'mViewSource', ...
        'Checked', 'on', 'Callback', 'kazan(''DataSource_Callback'',gcbo,[],guidata(gcbo))');
    handles.mViewOutput = uimenu(handles.mView, 'Label','Output', 'Tag', 'mViewOutput', ...
        'Callback', 'kazan(''DataSource_Callback'',gcbo,[],guidata(gcbo))');
    handles.mViewSnO = uimenu(handles.mView, 'Label','Source and Output', 'Tag', 'mViewSnO', ...
        'Callback', 'kazan(''DataSource_Callback'',gcbo,[],guidata(gcbo))');
    handles.mViewReal = uimenu(handles.mView, 'Label','Show real part', 'Tag', 'mViewReal', ...
        'Separator', 'on', 'Checked', 'on', 'Accelerator', 'R', ...
        'Callback', 'kazan(''DataType_Callback'',gcbo,[],guidata(gcbo))');
    handles.mViewImag = uimenu(handles.mView, 'Label','Show imaginary part', 'Tag', 'mViewImag', ...
        'Accelerator', 'I', 'Callback', 'kazan(''DataType_Callback'',gcbo,[],guidata(gcbo))');
    handles.mViewMagnitude = uimenu(handles.mView, 'Label','Show magnitude', 'Tag', 'mViewMagnitude', ...
        'Accelerator', 'M', 'Callback', 'kazan(''DataType_Callback'',gcbo,[],guidata(gcbo))');
    handles.mView1D = uimenu(handles.mView, 'Label','1D Slice', 'Tag', 'mView1D', ...
        'Separator', 'on', 'Callback', 'kazan(''PlotType_Callback'',gcbo,[],guidata(gcbo))');
    handles.mViewStacked = uimenu(handles.mView, 'Label','Stacked Plot', 'Tag', 'mViewStacked', ...
        'Checked', 'on', ...
        'Callback', 'kazan(''PlotType_Callback'',gcbo,[],guidata(gcbo))');
    handles.mViewContour = uimenu(handles.mView, 'Label','Contour Plot', 'Tag', 'mViewContour', ...
        'Callback', 'kazan(''PlotType_Callback'',gcbo,[],guidata(gcbo))');
    handles.mViewDensity = uimenu(handles.mView, 'Label','Density Plot', 'Tag', 'mViewDensity', ...
        'Callback', 'kazan(''PlotType_Callback'',gcbo,[],guidata(gcbo))');
    handles.mViewAxistight = uimenu(handles.mView, 'Label','Axis tight', 'Tag', 'mViewAxistight', ...
        'Checked', 'on', 'Separator', 'on', 'Accelerator', 'T',...
        'Callback', 'kazan(''mViewAxisSet'',gcbo,[],guidata(gcbo))');
    handles.mViewAxisauto = uimenu(handles.mView, 'Label','Axis auto', 'Tag', 'mViewAxisauto', ...
        'Accelerator', 'A', ...
        'Callback', 'kazan(''mViewAxisSet'',gcbo,[],guidata(gcbo))');
    handles.mViewAxisfixed = uimenu(handles.mView, 'Label','Axis fixed', 'Tag', 'mViewAxisfixed', ...
        'Accelerator', 'F', ...
        'Callback', 'kazan(''mViewAxisSet'',gcbo,[],guidata(gcbo))');
    handles.mViewShowGrid = uimenu(handles.mView, 'Label','Show Grid', 'Tag', 'mViewShowGrid', ...
        'Separator', 'on', 'Callback', 'kazan(''mViewAxisSet'',gcbo,[],guidata(gcbo))');
    handles.mViewAxisnormal = uimenu(handles.mView, 'Label','Aspect normaL', 'Tag', 'mViewAxisnormal', ...
        'Separator', 'on',  'Accelerator', 'L', ...
        'Callback', 'kazan(''mViewAxisSet'',gcbo,[],guidata(gcbo))');
    handles.mViewAxissquare = uimenu(handles.mView, 'Label','Aspect square', 'Tag', 'mViewAxissquare', ...
        'Accelerator', 'S', ...
        'Callback', 'kazan(''mViewAxisSet'',gcbo,[],guidata(gcbo))');
    handles.mViewAxisEqual = uimenu(handles.mView, 'Label','Axis equal', 'Tag', 'mViewAxisEqual', ...
        'Accelerator', 'Q', ...
        'Callback', 'kazan(''mViewAxisSet'',gcbo,[],guidata(gcbo))');
    handles.mViewRange = uimenu(handles.mView, 'Label','Select View Area', 'Tag', 'mViewRange', ...
        'Separator', 'on', 'Accelerator', 'V', ... 
        'Callback', 'kazan(''mViewRange'',gcbo,[],guidata(gcbo))');
    % menu Tools
    handles.mTools = uimenu(handles.MainFigure, 'Label','Tools', 'Tag', 'mTools');
    handles.mToolsClearSource = uimenu(handles.mTools, 'Label','Clear Source', 'Tag', 'mToolsClearSource', ...
        'Callback', 'kazan(''mTools'',gcbo,111,guidata(gcbo))');
    handles.mToolsFixProcessing = uimenu(handles.mTools, 'Label','Fix Source Processing', 'Tag', 'mToolsFixProcessing', ...
        'Callback', 'kazan(''mTools'',gcbo,111,guidata(gcbo))');
    handles.mToolsFixSlice = uimenu(handles.mTools, 'Label','Cut Slice', 'Tag', 'mToolsFixSlice', ...
        'Callback', 'kazan(''mTools'',gcbo,111,guidata(gcbo))');
    handles.mToolsLoadscript = uimenu(handles.mTools, 'Label','Show Data Load script', 'Tag', 'mToolsLoadscript', ...
        'Separator', 'on', 'Callback', 'kazan(''mTools'',gcbo,111,guidata(gcbo))');
    handles.mToolsExecLoadscript = uimenu(handles.mTools, 'Label','Execute Data Load Script', 'Tag', 'mToolsExecLoadscript', ...
        'Callback', 'kazan(''mTools'',gcbo,111,guidata(gcbo))');
    handles.mToolsCopyLoadscript = uimenu(handles.mTools, 'Label','Copy Load Script to Clipboard', 'Tag', 'mToolsCopyLoadscript', ...
        'Callback', 'kazan(''mTools'',gcbo,111,guidata(gcbo))');
    handles.mToolsShowinFigure = uimenu(handles.mTools, 'Label','Show in Separate Figure', 'Tag', 'mToolsShowinFigure', ...
        'Separator', 'on', 'Callback', 'kazan(''mTools'',gcbo,111,guidata(gcbo))');
    handles.mToolsAssociatedwindow = uimenu(handles.mTools, 'Label','No Figure Is Associated', 'Tag', 'mToolsAssociatedwindow');
    handles.mToolsAssociatedWindowReset = uimenu(handles.mToolsAssociatedwindow, 'Label','Reset', 'Tag', 'mToolsAssociatedWindowReset', ...
        'Callback', 'kazan(''mTools'',gcbo,[],guidata(gcbo))');
    handles.mToolsAssociatedWindowSet = uimenu(handles.mToolsAssociatedwindow, 'Label','To Set use KV_Tools/Associate of any KV-generated figure', 'Tag', 'mToolsAssociatedWindowSet', ...
        'Callback', '');
    handles.mToolsShowinwindow = uimenu(handles.mTools, 'Label','Show in Axis', 'Tag', 'mToolsShowinwindow');
    handles.mToolsShowinwindow0 = uimenu(handles.mToolsShowinwindow, 'Label','Open figure', 'Tag', 'mToolsShowinwindow0', ...
        'Callback', 'kazan(''mTools'',gcbo,[],guidata(gcbo))');
    handles.mToolsShowinwindow2 = uimenu(handles.mToolsShowinwindow, 'Label','1 / [1 2]', 'Tag', 'mToolsShowinwindow2', ...
        'Callback', 'kazan(''mTools'',gcbo,112,guidata(gcbo))');
    handles.mToolsShowinwindow3 = uimenu(handles.mToolsShowinwindow, 'Label','2 / [1 2]', 'Tag', 'mToolsShowinwindow3', ...
        'Callback', 'kazan(''mTools'',gcbo,212,guidata(gcbo))');
    handles.mToolsShowinwindow4 = uimenu(handles.mToolsShowinwindow, 'Label','1 / [2 1]', 'Tag', 'mToolsShowinwindow4', ...
        'Callback', 'kazan(''mTools'',gcbo,121,guidata(gcbo))');
    handles.mToolsShowinwindow5 = uimenu(handles.mToolsShowinwindow, 'Label','2 / [2 1]', 'Tag', 'mToolsShowinwindow5', ...
        'Callback', 'kazan(''mTools'',gcbo,221,guidata(gcbo))');
    handles.mToolsSkyline = uimenu(handles.mTools, 'Label','Show 2D density skyline plot', 'Tag', 'mToolsSkyline', ...
        'Callback', 'kazan(''mTools'',gcbo,[],guidata(gcbo))');
    handles.mToolsSkyline1 = uimenu(handles.mTools, 'Label','Show 2D contour skyline plot', 'Tag', 'mToolsSkyline1', ...
        'Callback', 'kazan(''mTools'',gcbo,[],guidata(gcbo))');
    handles.mToolsZoom = uimenu(handles.mTools, 'Label','Zoom', 'Tag', 'mToolsZoom', 'Separator', 'on', ...
        'Callback', 'kazan(''mTools'',gcbo,[],guidata(gcbo))', ...
        'Accelerator', 'Z');
    handles.mToolsAutoZoomOff = uimenu(handles.mTools, 'Label','One Zoom per Time', 'Tag', 'mToolsAutoZoomOff', ...
        'Callback', 'kazan(''mTools'',gcbo,[],guidata(gcbo))', ...
        'Checked', handles.autoZoomOff);
    handles.mToolsZoomFixed = uimenu(handles.mTools, 'Label','Set Fixed Axis', 'Tag', 'mToolsZoomFixed', ...
        'Callback', 'kazan(''mTools'',gcbo,[],guidata(gcbo))', ...
        'Checked', handles.FixedZoomOpt);
    handles.mCoord = uimenu(handles.mTools, 'Label','Measure', 'Tag', 'mCoord', 'Separator', 'on', ...
        'Callback', 'kazan(''mTools'',gcbo,[],guidata(gcbo))', ...
        'Accelerator', 'E');
    % menu Plugins
    handles.mPlugins = uimenu(handles.MainFigure, 'Label','Plugins', 'Tag', 'mPlugins', 'Separator', 'off');
    % menu Help
    handles.mHelp = uimenu(handles.MainFigure, 'Label','Help', 'Tag', 'mHelp', 'Separator', 'off');
    handles.mHelp3 = uimenu(handles.mHelp, 'Label','How to Start', 'Tag', 'mHelp', 'Separator', 'off', 'Callback', ...
        'kazan(''mHelp'',gcbo,[],guidata(gcbo))');
    handles.mHelp4 = uimenu(handles.mHelp, 'Label','Data Load Routines', 'Tag', 'mHelp', 'Separator', 'off', 'Callback', ...
        'kazan(''mHelp'',gcbo,[],guidata(gcbo))');
    handles.mHelp1 = uimenu(handles.mHelp, 'Label','Axis Data Structure (ax)', 'Tag', 'mHelp', 'Separator', 'off', 'Callback', ...
        'kazan(''mHelp'',gcbo,[],guidata(gcbo))');
    handles.mHelp2 = uimenu(handles.mHelp, 'Label','Mouse Cursor', 'Tag', 'mHelp', 'Separator', 'off', 'Callback', ...
        'kazan(''mHelp'',gcbo,[],guidata(gcbo))');
    handles.mHelp5 = uimenu(handles.mHelp, 'Label','Plugins', 'Tag', 'mHelp', 'Separator', 'off', 'Callback', ...
        'kazan(''mHelp'',gcbo,[],guidata(gcbo))');
    handles.mHelp6 = uimenu(handles.mHelp, 'Label','How to Print', 'Tag', 'mHelp', 'Separator', 'off', 'Callback', ...
        'kazan(''mHelp'',gcbo,[],guidata(gcbo))');
    handles.mHelpAbout = uimenu(handles.mHelp, 'Label','About Kazan viewer', 'Tag', 'mAbout', 'Separator', 'on', 'Callback', ...
        'kazan(''mAbout'',gcbo,[],guidata(gcbo))');
    %**********************************************************************
    %**********************************************************************
    %**********************************************************************

    %**********************************************************************
    %*****************   C o n t e x t     m e n u    *********************
    %**********************************************************************
    handles.cmContextMenu1 = uicontextmenu('Parent', handles.MainFigure);
    handles.cmNormal = uimenu(handles.cmContextMenu1, 'Label','Normal', 'Tag', 'cmNormal', ...
        'Checked', 'on', ...
        'Callback','kazan(''cmContexMenu1'',gcbo,[],guidata(gcbo))');
    handles.cmDiff0 = uimenu(handles.cmContextMenu1, 'Label','diff', 'Tag', 'cmDiff0', ...
        'Callback', 'kazan(''cmContexMenu1'',gcbo,[],guidata(gcbo))');
    handles.cmDiff1 = uimenu(handles.cmContextMenu1, 'Label','Ps.mod., (ps_mod)', 'Tag', 'cmDiff1', ...
        'Callback', 'kazan(''cmContexMenu1'',gcbo,[],guidata(gcbo))');
    handles.cmIntegrate = uimenu(handles.cmContextMenu1, 'Label','Integrate', 'Tag', 'cmIntegrate', ...
        'Callback', 'kazan(''cmContexMenu1'',gcbo,[],guidata(gcbo))');
    handles.cmNofilter = uimenu(handles.cmContextMenu1, 'Label','Not Filtered', 'Tag', 'cmNofilter', ...
        'Separator', 'on', 'Checked', 'on', ...
        'Callback',  'kazan(''Filter_CallBack'',gcbo,[],guidata(gcbo))');
    handles.cmRCfilter = uimenu(handles.cmContextMenu1, 'Label','RC (tc)', 'Tag', 'cmRCfilter', ...
        'Callback',  'kazan(''Filter_CallBack'',gcbo,[],guidata(gcbo))');
    handles.cmMvgaver = uimenu(handles.cmContextMenu1, 'Label','M.Aver (tc)', 'Tag', 'cmMvgaver', ...
        'Callback',  'kazan(''Filter_CallBack'',gcbo,[],guidata(gcbo))');
    handles.cmSavGol = uimenu(handles.cmContextMenu1, 'Label','Sav-Gol (tc, pol)', 'Tag', 'cmSavGol', ...
        'Callback',  'kazan(''Filter_CallBack'',gcbo,[],guidata(gcbo))');
    handles.cmApplyTo = uimenu(handles.cmContextMenu1, 'Label','Apply to:', 'Tag', 'cmApplyTo');
    %**********************************************************************
    handles.cmFiltAuto = uimenu(handles.cmApplyTo, 'Label','Auto', 'Tag', 'cmFiltAuto', ...
        'Callback',  'kazan(''Filter_CallBack'',gcbo,[],guidata(gcbo))');
    handles.cmFilt1D = uimenu(handles.cmApplyTo, 'Label','First Dim', 'Tag', 'cmFilt1D', ...
        'Callback',  'kazan(''Filter_CallBack'',gcbo,[],guidata(gcbo))');
    handles.cmFilt2D = uimenu(handles.cmApplyTo, 'Label','Second Dim', 'Tag', 'cmFilt2D', ...
        'Callback',  'kazan(''Filter_CallBack'',gcbo,[],guidata(gcbo))');
    handles.cmFilt1D2D = uimenu(handles.cmApplyTo, 'Label','Both', 'Tag', 'cmFilt1D2D', ...
        'Callback',  'kazan(''Filter_CallBack'',gcbo,[],guidata(gcbo))');
    %**********************************************************************
    handles.cmxy = uimenu(handles.cmContextMenu1, 'Label','x-y', 'Tag', 'cmxy', ...
        'Separator', 'on', 'Callback', ...
        'kazan(''cmContexMenu2'',gcbo,[],guidata(gcbo))');
    handles.xxxxxx = uimenu(handles.cmContextMenu1, 'Label','X', 'Tag', 'xxxxxx');
    %**********************************************************************
    handles.cmxLinear = uimenu(handles.xxxxxx, 'Label', 'X', 'Tag', 'cmxLinear', ...
        'Checked', 'on', 'Callback', 'kazan(''cmContexMenu2'',gcbo,[],guidata(gcbo))');
    handles.cmxPoints = uimenu(handles.xxxxxx, 'Label', 'X Points', 'Tag', 'cmxPoints', ...
        'Callback', 'kazan(''cmContexMenu2'',gcbo,[],guidata(gcbo))');
    handles.cmReverseX = uimenu(handles.xxxxxx, 'Label', 'Reverse X', 'Tag', 'cmReverseX', ...
        'Callback', 'kazan(''cmContexMenu2'',gcbo,[],guidata(gcbo))');
    handles.cmSemilogX = uimenu(handles.xxxxxx, 'Label', 'X Log', 'Tag', 'cmSemilogX', ...
        'Callback', 'kazan(''cmContexMenu2'',gcbo,[],guidata(gcbo))');
    handles.cmGspace = uimenu(handles.xxxxxx, 'Label','X g-space (freq1)', 'Tag', 'cmGspace', ...
        'Callback', 'kazan(''cmContexMenu2'',gcbo,[],guidata(gcbo))');
    handles.cmxPPM = uimenu(handles.xxxxxx, 'Label','X PPM (reference)', 'Tag', 'cmxPPM', ...
        'Callback', 'kazan(''cmContexMenu2'',gcbo,[],guidata(gcbo))');
    handles.cmxPPMmin1 = uimenu(handles.xxxxxx, 'Label','X PPM-1 (reference)', 'Tag', 'cmxPPMmin1', ...
        'Callback', 'kazan(''cmContexMenu2'',gcbo,[],guidata(gcbo))');
    %**********************************************************************
    handles.yyyyyy = uimenu(handles.cmContextMenu1, 'Label','Y', 'Tag', 'yyyyyy');
    %**********************************************************************
    handles.cmyLinear = uimenu(handles.yyyyyy, 'Label', 'Y', 'Tag', 'cmyLinear', ...
        'Checked', 'on', 'Callback', 'kazan(''cmContexMenu2'',gcbo,[],guidata(gcbo))');
    handles.cmyPoints = uimenu(handles.yyyyyy, 'Label', 'Y Points', 'Tag', 'cmyPoints', ...
        'Callback', 'kazan(''cmContexMenu2'',gcbo,[],guidata(gcbo))');
    handles.cmReverseY = uimenu(handles.yyyyyy, 'Label', 'Reverse Y', 'Tag', 'cmReverseY', ...
        'Callback', 'kazan(''cmContexMenu2'',gcbo,[],guidata(gcbo))');
    handles.cmSemilogY = uimenu(handles.yyyyyy, 'Label', 'Y Log', 'Tag', 'cmSemilogY', ...
        'Callback', 'kazan(''cmContexMenu2'',gcbo,[],guidata(gcbo))');
    handles.cmyPPM = uimenu(handles.yyyyyy, 'Label','Y PPM (reference)', 'Tag', 'cmyPPM', ...
        'Callback', 'kazan(''cmContexMenu2'',gcbo,[],guidata(gcbo))');
    handles.cmyPPMmin1 = uimenu(handles.yyyyyy, 'Label','Y PPM-1 (reference)', 'Tag', 'cmyPPMmin1', ...
        'Callback', 'kazan(''cmContexMenu2'',gcbo,[],guidata(gcbo))');
    %**********************************************************************
    %**********************************************************************
    %**********************************************************************
    
    % left panel
    handles.fDirBrowser = uicontrol(handles.MainFigure,'Style', 'frame', 'Tag', 'fDirBrowser', ...
        'Units','characters');
    handles.eDirectory = uicontrol(handles.MainFigure,'Style', 'edit', 'Tag', 'eDirectory', ...
        'Units','characters', 'Callback', 'kazan(''eDirectory_Callback'',gcbo,[],guidata(gcbo))');
    handles.pmDatatype = uicontrol(handles.MainFigure,'Style', 'popupmenu', 'Tag', 'pmDatatype', ...
        'Units','characters', 'String', '*.*', 'Callback', 'kazan(''pmDatatype_Callback'',gcbo,[],guidata(gcbo))');
    handles.pbUp = uicontrol(handles.MainFigure,'Style', 'pushbutton', 'Tag', 'pbUp', ...
        'Units','characters', 'String', 'Up', 'Callback', 'kazan(''pbUp_Callback'',gcbo,[],guidata(gcbo))');
    handles.lDirlist = uicontrol(handles.MainFigure,'Style', 'listbox', 'Tag', 'lDirlist', ...
        'Units','characters', 'String', 'ss', 'Callback', 'kazan(''lDirlist_Callback'',gcbo,[],guidata(gcbo))');
    handles.pmExt = uicontrol(handles.MainFigure,'Style', 'popupmenu', 'Tag', 'pmExt', ...
        'Units','characters', 'String', '*.*', 'Callback', 'kazan(''pmExt_Callback'',gcbo,[],guidata(gcbo))');

    set(handles.eDirectory, 'TooltipString', 'dir');

    % Upper panel
    handles.fGrOptions = uicontrol(handles.MainFigure,'Style', 'frame', 'Tag', 'fGrOptions', 'Units','characters');
    handles.rbSource = uicontrol(handles.MainFigure,'Style', 'toggle', 'Tag', 'rbSource', 'Value', 1, 'Units','characters', ...
        'String', 'Source', 'Callback', 'kazan(''SetDataSource'',1,guidata(gcbo))');
    handles.rbOutput = uicontrol(handles.MainFigure,'Style', 'toggle', 'Tag', 'rbOutput', 'Value', 0, 'Units','characters',...
        'String', 'Output', 'Callback', 'kazan(''SetDataSource'',2,guidata(gcbo))');
    handles.rbOutToSrc = uicontrol(handles.MainFigure,'Style', 'pushbutton', 'Tag', 'rbOutToSrc', 'Units','characters',...
        'String', ' << ', 'Callback', 'kazan(''mFile2Source'',gcbo,[],guidata(gcbo))');

    handles.rbReal = uicontrol(handles.MainFigure,'Style', 'toggle', 'Tag', 'rbReal', 'Units','characters', ...
        'String', 'Real', 'Callback', 'kazan(''SetDataPart'',1,guidata(gcbo))', 'Value', 1);
    handles.rbImag = uicontrol(handles.MainFigure,'Style', 'toggle', 'Tag', 'rbImag', 'Units','characters', ...
        'String', 'Re/Im', 'Callback', 'kazan(''SetDataPart'',2,guidata(gcbo))', 'Value', 0);
    handles.rbMagnitude = uicontrol(handles.MainFigure,'Style', 'toggle', 'Tag', 'rbMagnitude', 'Units','characters', ...
        'String', 'Magn', 'Callback', 'kazan(''SetDataPart'',3,guidata(gcbo))', 'Value', 0);

    % plot type buttons
    handles.tb1D = uicontrol(handles.MainFigure,'Style', 'toggle', 'Tag', 'tb1D', 'Units','characters', ...
        'String', '1D', 'Callback', 'kazan(''SetPlotType'',1,guidata(gcbo))', 'Value', 1);
    handles.tbStacked = uicontrol(handles.MainFigure,'Style', 'toggle', 'Tag', 'tbStacked', 'Units','characters', ...
        'String', 'Stack', 'Callback', 'kazan(''SetPlotType'',2,guidata(gcbo))', 'Value', 0);
    handles.tbContour = uicontrol(handles.MainFigure,'Style', 'toggle', 'Tag', 'tbContour', 'Units','characters', ...
        'String', 'Contour', 'Callback', 'kazan(''SetPlotType'',3,guidata(gcbo))', 'Value', 0);
    handles.tbDensity = uicontrol(handles.MainFigure,'Style', 'toggle', 'Tag', 'tbDensity', 'Units','characters', ...
        'String', 'Density', 'Callback', 'kazan(''SetPlotType'',4,guidata(gcbo))', 'Value', 0);

    % zoom buttons     
    handles.tbMeasure = uicontrol(handles.MainFigure,'Style', 'pushbutton', 'Tag', 'tbMeasure', 'Units','characters', ...
        'String', '.', 'Callback', 'kazan(''mTools'',gcbo,[],guidata(gcbo))', 'Value', 0, 'ToolTip', 'Measure');
    handles.tbX_half = uicontrol(handles.MainFigure,'Style', 'pushbutton', 'Tag', 'tbX_half', 'Units','characters', ...
        'String', '><', 'Callback', 'kazan(''Zoom_Callback'',gcbo,[],guidata(gcbo))', 'Value', 0, 'ToolTip', 'Zoom out X');
    handles.tbX_twice = uicontrol(handles.MainFigure,'Style', 'pushbutton', 'Tag', 'tbX_twice', 'Units','characters', ...
        'String', '<>', 'Callback', 'kazan(''Zoom_Callback'',gcbo,[],guidata(gcbo))', 'Value', 0, 'ToolTip', 'Zoom in X');
    handles.tbY_half = uicontrol(handles.MainFigure,'Style', 'pushbutton', 'Tag', 'tbX_half', 'Units','characters', ...
        'String', 'X', 'Callback', 'kazan(''Zoom_Callback'',gcbo,[],guidata(gcbo))', 'Value', 0, 'ToolTip', 'Zoom out Y');
    handles.tbY_twice = uicontrol(handles.MainFigure,'Style', 'pushbutton', 'Tag', 'tbX_twice', 'Units','characters', ...
        'String', '^', 'Callback', 'kazan(''Zoom_Callback'',gcbo,[],guidata(gcbo))', 'Value', 0, 'ToolTip', 'Zoom in Y');
    handles.tbImReX = uicontrol(handles.MainFigure,'Style', 'pushbutton', 'Tag', 'tbImReX', 'Units','characters', ...
        'String', 'xRI', 'Callback', 'kazan(''Zoom_Callback'',gcbo,[],guidata(gcbo))', 'Value', 0, 'ToolTip', 'xRe=xIm');    
%%%%%%%%%%%%%%%%%%%%%%% for ML >7
    if MatlabVer >= 7,
        handles.tbPan = uicontrol(handles.MainFigure,'Style', 'toggle', 'Tag', 'tbPan', 'Units','characters', ...
        'String', 'W', 'Callback', 'kazan(''Zoom_Callback'',gcbo,[],guidata(gcbo))', 'Value', 0, 'ToolTip', 'Pan');  
    else
        handles.tbPan = 0;
    end
    
    % Main part
    handles.tUnits = uicontrol(handles.MainFigure,'Style', 'text', 'Tag', 'tUnits', 'Units','characters');
    handles.stCoordinates = uicontrol(handles.MainFigure,'Style', 'text', 'Tag', 'stCoordinates', 'Units','characters');
    handles.aReal = axes('Parent', handles.MainFigure, 'Tag', 'aReal', 'Units','characters', 'UIContextMenu', handles.cmContextMenu1);
    handles.aImag = axes('Parent', handles.MainFigure, 'Tag', 'aImag', 'Units','characters', 'UIContextMenu', handles.cmContextMenu1);
    handles.pmSelPar = uicontrol(handles.MainFigure,'Style', 'popupmenu', 'Tag', 'pmSelPar', 'Units','characters', ...
            'Callback', 'kazan(''smSelPar_Callback'',gcbo,[],guidata(gcbo))');
    handles.ePar =  uicontrol(handles.MainFigure,'Style', 'edit', 'Tag', 'ePar', 'Units','characters', ...
            'Callback', 'kazan(''ePar_Callback'',gcbo,[],guidata(gcbo))');
    handles.sl2Dplot  = uicontrol(handles.MainFigure,'Style', 'slider', 'Tag', 'sl2Dplot', 'Visible', 'off', ...
            'Units','characters', 'Callback', 'kazan(''sl2Dplot_Callback'',gcbo,[],guidata(gcbo))');
    handles.pb2DSelect  = uicontrol(handles.MainFigure,'Style', 'pushbutton', 'Tag', 'pb2DSelect', 'Visible', 'off', 'String', 'X', ...
            'Units','characters', 'Callback', 'kazan(''pb2DSelect_Callback'',gcbo,[],guidata(gcbo))');

    handles.fMoveDataBr = uicontrol(handles.MainFigure,'Style', 'frame', 'Tag', 'fMoveDataBr', ...
        'Callback', 'kazan(''fMoveDataBr_Callback'',gcbo, guidata(gcbo))','Units','characters', ...
        'ButtonDownFcn', 'kazan(''MainFigure_ButtonDownFcn'', gcbo, [], guidata(gcbo))', 'String', '', ...
        'Enable', 'off');
    
    %   KAZAN VIEWER logo
    handles.src.ax.xlabel = 'KAZAN viewer logo,';
    handles.src.ax.type = 'data';
    handles.src.ax.filt = '';
    handles.src.ax.diff = '';
    sz = 0.3;
    kx = [-.15 0 0 0 1 0 1]*1.3; ky = [1 1 0 .5 1 .5 0];
    ax = [0 .5 .75 .25 .75 1]; ay = [0 1 .5 .5 .5 0];
    zx = [0 1 0  0  0 1 0 1 1 1]; zy = [0 1 1 .92 1 1  0 0 .08 0];
    nx = [0 0 1 1 1+.15]*1.4; ny = [0 1 0  1  1];
    handles.src.ax.x = [8,8,-2,-2, kx-.3,ax+1.2,zx+2.4,ax+3.6,nx+4.8, ...
        8,8,6.2,8,8,-2,-2,-0.3,-2,-2]'*sz;
    handles.src.ax.y = 1;
    handles.src.y = [1,2,2,1,ky,ay,zy,ay,ny, ...
        1,0,0,0,-1,-1,0,0,0,1]';
    handles.lscript = [];
    handles.DataSource = 1;
    handles.DataPart = 1;
    handles.PlotType = 1;
    handles.xrepresent = 'plot';
    handles.yrepresent = 'plot';
    handles.axisset = 'tight';
    handles.axisaspect = 'normal';
    handles.axisaspect1 = 'normal';
    handles.axisgrid = 'off';
    %   The test dataset
    x = (1:200)';
    handles.out{1}.xlabel = 'Frequency, MHz';
    handles.out{1}.ax.x =  x;
    handles.out{1}.ax.y =  1;
    handles.out{1}.ax.type = 'data';
    handles.out{1}.y = exp(-x/26).*(sin(x*.3)+1i*cos(x*.4)*.5);
    %   Viewer options
    handles.PluginHandle = [];
    handles.ReadOptHandle = [];
    handles.ReadOptLast = '';
    handles.SecondCopy = 1;
    handles.ext = '';
    handles.Datatype = '';
    handles.fname = 'none';
    handles.fhistory = {};
    handles.dhistory = {};
    handles.favordir = {};
    handles.process = '';
    handles.filter = '';
    handles.isScroller = 0;
    handles.Projection = 1;
    handles.AntiProjection = 2;
    handles.Selection  = 1;
    handles.DataPartSize = 40;
    handles.AssociatedFigure = [];
    handles.changePart = 0;
    str = {'title','xlabel','ylabel','Color','LineStyle','LineWidth','Marker','dx','s','freq1', 'reference', 'cf', 'ps_mod', 'ps_mod_harm', 'tc', 'pol', 'contour','yslicelabel'};
    set(handles.pmSelPar, 'String', str, 'Value', 1)
    handles.prop.values = {'?','?','?','b','-','1','none','0', '0', '0', '0', '0', '5', '1','1', '2', '.2:.1:.6',''};
    handles.prop.types = {'s','s','s','s','s','f','s','f', 'f', 'f', 'f', 'f', 'f', 'f','f', 'f', 'f','s'};


    handles.DirStartCode =' ['; %<HTML><FONT color=#763F18"><b>
    handles.DirEndCode=']';
    
    % Data extensions
    ext_list{1} = '*.*';
    ext_list{end+1} = '*.aqs';
    ext_list{end+1} = '*.bkn';
    ext_list{end+1} = '*.jpg';
    ext_list{end+1} = '*.dat';
    ext_list{end+1} = '*.d00';
    ext_list{end+1} = '*.d01';
    ext_list{end+1} = '*.dsc';
    ext_list{end+1} = '*.exp';
    ext_list{end+1} = '*.mos';
    ext_list{end+1} = '*.par';
    ext_list{end+1} = '*.rmn';
    ext_list{end+1} = '*.txt';
    ext_list{end+1} = '*.scr';
    ext_list{end+1} = '*.sp';
    ext_list{end+1} = '*.wcw';
    ext_list{end+1} = '*.w2p';    
    ext_list{end+1} = '*.wtr';
    ext_list{end+1} = '*.xyz';
    ext_list{end+1} = '*.res';
    ext_list{end+1} = '*.xml';
    ext_list{end+1} = '*.spa';
    ext_list{end+1} = '*.cmbl';
    ext_list{end+1} = '*.';
    set(handles.pmExt, 'String',ext_list);
   


    set(handles.MainFigure, 'ResizeFcn', 'kazan(''MainFigure_ResizeFcn'', gcbo,[],guidata(gcbo))');
    set(handles.MainFigure, 'Resize', 'on');

    % load directory
    set(handles.eDirectory, 'String', handles.dir_path, 'TooltipString', handles.dir_path)
    %     set(handles.eDirectory, 'TooltipString', 'dir');
    handles = pmExt_Callback(gcbo, [], handles, varargin);

    set(handles.MainFigure, 'CreateFcn', 'kazan(''MainFigure_CreateFcn'', gcbo,[],guidata(gcbo))');
    set(handles.MainFigure, 'DeleteFcn', 'kazan(''MainFigure_DeleteFcn'', gcbo,[],guidata(gcbo))');
    set(handles.MainFigure, 'WindowButtonDownFcn', 'kazan(''MainFigure_ButtonDownFcn'', gcbo,[],guidata(gcbo))');
    set(handles.MainFigure, 'WindowButtonUpFcn', 'kazan(''MainFigure_ButtonUpFcn'', gcbo,[],guidata(gcbo))');
    set(handles.MainFigure, 'WindowButtonMotionFcn', 'kazan(''MainFigure_ButtonMotionFcn'', gcbo,[],guidata(gcbo))');


    guidata(fig, handles);

    % load toolboxes for option 3
    if opc==3
        % alsi 14.03.2012 .... get plugins to a foulder
        mPluginsReload(handles.MainFigure, [], handles);
        handles = guidata(handles.MainFigure);
    end

    if nargout > 0
        varargout{1} = fig;
    end

    if isempty(oldfig),
        MainColor =  [0.8314, 0.8157, 0.7843];
        set(handles.MainFigure, 'Name', 'KAZAN viewer');
    else
        MainColor = [0.76, 0.734, 0.682];
        set(handles.MainFigure, 'Name',sprintf('[%d] KAZAN viewer',length(oldfig)+1));
    end
    set(handles.MainFigure, 'Color',MainColor);
    MainFigure_ResizeFcn(0, 0, handles);
    hdl = get(handles.MainFigure, 'Children');
    for k=1:length(hdl)
        if strcmp(get(hdl(k), 'Type'), 'uicontrol')
            set(hdl(k), 'BackgroundColor', MainColor);
        end
    end
    SetPlotType(2, handles);
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
%| 'MainFigure_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.MainFigure, handles.slider2. This
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

% --------------------------------------------------------------------
function varargout = MainFigure_CreateFcn(h, eventdata, handles, varargin)
varargout = {};

% --------------------------------------------------------------------
function mAbout(h, eventdata, handles)
msgbox({'KAZAN data viewer and plugins';'Version 2.0 beta (07.03.2013)'; 'by Boris Epel and Alexey Silakov, 2003-09'}, ...
    'About', 'help')

% --------------------------------------------------------------------
function fig_size = check_figure_size(fig)
% manual resize of the figure
un = get(0, 'Units');
set(0, 'Units', 'characters');
screen_size = get(0, 'ScreenSize');
set(0, 'Units', un);
fig_size = get(fig, 'Position');

% minimum size check
update = 0;
if fig_size(3) < 100, fig_size(3) = 100; update = 1; end
if fig_size(4) < 20, fig_size(4) = 20;   update = 1; end
if fig_size(2) + fig_size(4) > screen_size(4)
    fig_size(2) = screen_size(4) - fig_size(4); update = 1;
end

if update
    %% set(handles.MainFigure, 'Position', fig_size);
end
    
% --------------------------------------------------------------------
function handles = CreateFavorites(handles)
menus = get(handles.mFavorites, 'Children');
delete(findobj(menus, 'UserData', 'favorites'));
outstr = inimanaging(handles.inifile);

favorDir = {};
try  favorDir = outstr.KazanViewer.favorDir; catch end
n = 1;
callback = get(handles.mFavorAdd, 'Callback');
check = {'off' 'on'};
if ~isempty(favorDir)
    for ci = 1:length(favorDir)
        sep = check{(n == 1) + 1};
        uimenu(handles.mFavorites, 'Label',favorDir{ci}, 'Callback', callback, ...
            'Tag', ['mFavorites' n], 'UserData', 'favorites', 'Separator', sep);
        n = n + 1;
    end
end
% --------------------------------------------------------------------

function MainFigure_DeleteFcn(h, eventdata, handles, varargin)
delete(handles.PluginHandle);
try
    outstr = inimanaging(handles.inifile);
    outstr.KazanViewer.lastDir = get(handles.eDirectory, 'String');
    outstr.KazanViewer.autoZoomOff = handles.autoZoomOff;
    outstr.KazanViewer.FixedZoomOpt = handles.FixedZoomOpt;
    inimanaging(handles.inifile, outstr);
catch
end
% --------------------------------------------------------------------

function MainFigure_ResizeFcn(h, eventdata, handles, varargin)

if ~isfield(handles, 'aReal')
    return;
end

% Obtain size of the window
rectW = check_figure_size(handles.MainFigure);

% rectW - window size
%   DataSize                     GraphSize                   
% -----------------------------------------------------------
% |           ||               button_panel1                 | 
% |           ||---------------------------------------------|
% |           ||               button_panel2                 |
% |           ||---------------------------------------------|
% |           ||                                             |
% |           ||                                             |
% |           ||                axis_panel                   |
% |           ||                                             |
% |           ||                                             |
% |           ||                                             |
% |           ||---------------------------------------------|
% |           ||               parameters_panel              |
% -----------------------------------------------------------


Border = 0.4;
ButtonPanel_Height = 2.5;
Slider_Width = 1;
DataSize = safeget(handles,'DataPartSize',40); % in characters
GraphSize = rectW(3) - DataSize - Slider_Width;

data_panel     = [0,0,DataSize,rectW(4)];
slider         = [data_panel(1)+data_panel(3), 0,Slider_Width,rectW(4)];
graphics_panel = [slider(1)+slider(3),         0,GraphSize,rectW(4)];
button_panel1  = [graphics_panel(1), rectW(4)-ButtonPanel_Height, GraphSize, ButtonPanel_Height];
button_panel2  = [graphics_panel(1), rectW(4)-1.75*ButtonPanel_Height, GraphSize, ButtonPanel_Height];
parameters_panel  = [graphics_panel(1), 0, GraphSize, ButtonPanel_Height];
axis_panel     = graphics_panel + [0,ButtonPanel_Height,0,-2.7*ButtonPanel_Height];

% update of Graphics panel
rbSize = 9.1;
set(handles.fGrOptions, 'Position', button_panel1+Border*[1,1,-2,-2]);
rect_button = button_panel1 + Border*[2,1.25,0,-2.5];
rect_button(3) = rbSize-Border;
set(handles.rbSource, 'Position', rect_button); % button Source
rect_button(1) = rect_button(1) + rbSize;
rect_button(3) = rbSize*0.55-Border;
set(handles.rbOutToSrc, 'Position', rect_button); % button Output to Source
rect_button(1) = rect_button(1) + rbSize*0.55;
rect_button(3) = rbSize-Border;
set(handles.rbOutput, 'Position', rect_button); % button Output
rect_button(1) = rect_button(1) + rbSize*1.2;
set(handles.rbReal, 'Position', rect_button);
rect_button(1) = rect_button(1) + rbSize;
set(handles.rbImag, 'Position', rect_button);
rect_button(1) = rect_button(1) + rbSize;
set(handles.rbMagnitude, 'Position', rect_button);
rect_button(1) = rect_button(1) + rbSize*1.2;
set(handles.tb1D, 'Position', rect_button);
rect_button(1) = rect_button(1) + rbSize;
set(handles.tbStacked, 'Position', rect_button);
rect_button(1) = rect_button(1) + rbSize;
set(handles.tbContour, 'Position', rect_button);
rect_button(1) = rect_button(1) + rbSize;
set(handles.tbDensity, 'Position', rect_button);

rbSize = 5.0;
ButtonSpace = Border - .25;
rect_button = button_panel2 + Border*[2,1.25,0,-2.5];
rect_button(1) = button_panel2(1)+button_panel2(3)-rbSize;
rect_button(3) = rbSize-ButtonSpace;
set(handles.tbY_half, 'Position', rect_button);
rect_button(1) = rect_button(1) - rbSize;
set(handles.tbY_twice, 'Position', rect_button);
rect_button(1) = rect_button(1) - rbSize;
set(handles.tbX_half, 'Position', rect_button);
rect_button(1) = rect_button(1) - rbSize;
set(handles.tbX_twice, 'Position', rect_button);
rect_button(1) = rect_button(1) - rbSize;
set(handles.tbImReX, 'Position', rect_button);
rect_button(1) = rect_button(1) - rbSize;
set(handles.tbMeasure, 'Position', rect_button);
% if handles.tbPan~=0 % ONLY FOR MATLAB 7 OR HIGHER
    rect_button(1) = rect_button(1) - rbSize;
    set(handles.tbPan, 'Position', rect_button);
% end

set(handles.stCoordinates,'Position', button_panel2+[Border+14, Border - 0.2, -7*rbSize-14-Border, -2*Border-.2]);

% draw axis
isScroller = handles.isScroller;
Scroller_Width = 4;
ScrollerPlace = (1+Scroller_Width) * (isScroller > 0) + Border;
Border = Border * 2;
TextBorderX = Border + 7;
TextBorderY = 1.7;

switch handles.DataPart
    case {2}
        WindowHeight = axis_panel(4)/2;
        rectG = axis_panel + [TextBorderX, TextBorderY, -ScrollerPlace-TextBorderX-Border, -TextBorderY-WindowHeight];
        set(handles.aImag, 'Visible', 'on');
        set(handles.aImag, 'Position', rectG);
        rectG(2) = rectG(2) + WindowHeight;
        set(handles.aReal, 'Position', rectG);
    case {1,3}
        set(handles.aImag, 'Visible', 'off');
        rectG = axis_panel + [TextBorderX, TextBorderY, -ScrollerPlace-TextBorderX-Border, -TextBorderY];
        set(handles.aReal, 'Position', rectG);
end

% parameters position
rect_button = parameters_panel + [Border,.5*Border,0,-Border];
rect_button(3) = 15 - Border;
set(handles.pmSelPar, 'Position', rect_button);
rect_button(1) = rect_button(1) + 15;
rect_button(3) = 35 - Border;
set(handles.ePar, 'Position', rect_button);
% unit position
rect_button(1) = rect_button(1) + 35;
rect_button(3) = parameters_panel(1)+parameters_panel(3) - rect_button(1) - ScrollerPlace - 7;
set(handles.tUnits, 'Position', rect_button);
handles = DataLoadOptions(handles);
if handles.isScroller,
    rectSl([2,4]) = axis_panel([2,4]) + 1 * [1,-1]; 
    rectSl(1) = axis_panel(1) + axis_panel(3)- Scroller_Width - Border; 
    rectSl(3) = Scroller_Width;
    set(handles.sl2Dplot, 'Position', rectSl);
    rectSel = [rectSl(1), rectSl(2)-3, Scroller_Width, 2];
    set(handles.pb2DSelect, 'Position', rectSel);
end

% --------------------------------------------------------------------
function handles = pmExt_Callback(h, eventdata, handles, varargin)
str = get(handles.pmExt, 'String');
sel = get(handles.pmExt, 'Value');
if isempty(sel)
    if size(str,1) > 0, sel = 1;
    else return;
    end
end
switch str{sel}
    case '*.*'
        types{1}='ASCII';
        types{end+1}='ASCII2D';
        types{end+1}='UofC_bin';
        types{end+1}='RMN_1D';
        types{end+1}='RMN_2D';
        types{end+1}='TecMag';
        types{end+1}='XWINNMR';
        types{end+1}='WINNMR';
        types{end+1}='WEBMOSS';
        types{end+1}='MOSS';
        types{end+1}='OPUS_FTIR';
        types{end+1}='PerkinIR';
        types{end+1}='image';
        types{end+1}='SXYZ';
        types{end+1}='Wband_Brln';
        types{end+1}='ElChem_PAR';
        types{end+1}='SSRL-EXAFS';
        types{end+1}='CST';
        types{end+1}='AktaStart';
        types{end+1}='Agilent_BKN';
        types{end+1}='ESRxml';
        types{end+1}='NicoletSPA';
    case '*.'
        types{1}='ASCII';
        types{end+1}='ASCII2D';
        types{end+1}='RMN_1D';
        types{end+1}='RMN_2D';
        types{end+1}='Spinsight';
        types{end+1}='TecMag';
        types{end+1}='XWINNMR';
        types{end+1}='WINNMR';
    case '*.aqs'
        types{1}='WINNMR';
    case '*.dat'
        types{1}='ASCII';
        types{end+1}='ASCII2D';
        types{end+1}='XEPR';
        types{end+1}='RMNdat';
        types{end+1}='RMNdat2D';
        types{end+1}='transient';
        types{end+1}='CST';
    case '*.txt'
        types{1}='ASCII';
        types{end+1}='ASCII2D';
        types{end+1}='CST';
    case {'*.exp'}
        types{1}='WIS';
        types{2}='SpecMan';
    case {'*.par'}
        types{1}='XEPR';
        types{2}='XEPR_JSS';
        types{3}='ElChem_PAR';
    case {'*.dsc'}
        types{1}='XEPR';
        types{2}='XEPR_JSS';
    case '*.rmn'
        types{1}='RMN_1D';
        types{end+1}='RMN_2D';
    case {'*.d00', '*.exp'}
        types{1}='WIS';
    case '*.d01'
        types{1}='SpecMan';
    case '*.jpg'
        types{1}='image';
    case '*.xyz'
        types{1}='SXYZ';
    case '*.scr'
        types{1}='HYPERCH';  
    case {'*.wcw', '*.w2p', '*.wtr'}
        types{1}='Wband_Brln';
    case {'*.par'}
        types{1}='ElChem_PAR';    
    case {'*.mos'}
        types{1} = 'MOSS';
    case {'*.sp'}        
        types{1}='PerkinIR';
    case {'*.res'}
        types{1}='AktaStart';
    case {'*.bkn'}
        types{1} = 'Agilent_BKN';
    case {'*.xml'}
        types{1} = 'ESRxml';     
    case {'*.spa'}
        types{1} = 'NicoletSPA';     
    case {'*.cmbl'}
        types{1} = 'VernierCMBL';  
    otherwise
        types{1} = 'no types';
end
set(handles.pmDatatype, 'String', types, 'Value', 1);
handles.ext = str{sel};
handles = pmDatatype_Callback(h, eventdata, handles, varargin);
handles = LoadDirectory('', handles);
guidata(handles.MainFigure, handles);

% --------------------------------------------------------------------
% this function prepares correct place for different load options
function handles = DataLoadOptions(handles)
rectW = get(handles.MainFigure, 'Position');

Border = 0.4;
DataSize = handles.DataPartSize; % in characters

% update of Data read panel
rectW(4) = max(rectW(4), 20);
rectD = [Border, Border, max(DataSize-2*Border,0), max(rectW(4)-2*Border,0)];
set(handles.fDirBrowser, 'Position', rectD);
set(handles.fMoveDataBr, 'Position', [DataSize, Border, 1, rectW(4)-2*Border]);
controlSize = 1.6;
controlUpSize = 5.0;
InnerBorder = Border + .3;
rectDir = [rectD(1) + InnerBorder, rectD(2)+rectD(4)-Border-controlSize, rectD(3)-2*InnerBorder, controlSize];
set(handles.eDirectory, 'Position', rectDir);
rectDir(2) = rectDir(2) - controlSize-Border;
rectDir(3) = (DataSize-4*InnerBorder)/2;
set(handles.pmDatatype, 'Position', rectDir);
rectDir(1) = rectDir(1) + DataSize/2-InnerBorder;
rectDir(3) = rectDir(3)-controlUpSize;
set(handles.pmExt, 'Position', rectDir);
rectDir(1) = DataSize-InnerBorder-controlUpSize;
rectDir(3) = controlUpSize-InnerBorder;
set(handles.pbUp, 'Position', rectDir);

% DataSize = handles.DataPartSize; % in characters

InnerBorder = .3 + Border; controlHeight = 1.6;


rectD = [Border, Border, DataSize-2*Border, rectW(4)-2*Border];
set(handles.fDirBrowser, 'Position', rectD);
rectD = [Border + InnerBorder, Border+InnerBorder, ...
    DataSize-2*(Border+InnerBorder), ...
    rectW(4)-2*(Border+InnerBorder+controlHeight)-Border];
% set(handles.lDirlist, 'Position', rectD);

% Long control
controlLeft  = rectD(1);
controlWidth = DataSize-2*Border-2*InnerBorder;
controlTop = rectD(2) + (0:3)*(controlHeight + Border);

% Two controls
controlLeft2  = controlLeft + controlWidth/2 + InnerBorder;
controlWidth2 = controlWidth/2 - InnerBorder;

% One description and two controls
controlLeft32  = controlLeft + controlWidth/3 + InnerBorder;
controlLeft33  = controlLeft + 2*controlWidth/3 + InnerBorder;
controlWidth3  = (controlWidth-2*InnerBorder)/3 - InnerBorder;

otype = [handles.ext handles.Datatype];
if ~strcmp(handles.ReadOptLast, otype)
    delete(handles.ReadOptHandle)
    handles.ReadOptHandle = [];
end
handles.ReadOptLast = otype;
switch handles.ReadOptLast
    case {'*.datASCII', '*.*ASCII', '*.txtASCII', '*.ASCII'}
        spaceneeded = 2*controlHeight+2*Border;

        % create three controls for unit, delimiter and read format
        if isempty(handles.ReadOptHandle), % create
            handles.ReadOptHandle(1) = uicontrol(handles.MainFigure, 'Style', 'popupmenu',...
                'String', {'tab', 'space', ',', ';'}, ...
                'UserData', {'\t', ' ', ',', ';'}, 'Units', 'characters', ...
                'Position', [controlLeft,controlTop(1),controlWidth2,controlHeight], ...
                'Callback', 'kazan(''lDirlist_Callback'',gcbo,[],guidata(gcbo))');
            handles.ReadOptHandle(2) = uicontrol(handles.MainFigure, 'Style', 'popupmenu',...
                'String', {'Y', 'rY, iY', 'X,Y', 'X,rY,iY', 'idx, X, Y', 'idx, X, rY, iY'}, 'Units', 'characters', ...
                'Value', 4,'Position', [controlLeft2,controlTop(1),controlWidth2,controlHeight], ...
                'Callback', 'kazan(''lDirlist_Callback'',gcbo,[],guidata(gcbo))');
            handles.ReadOptHandle(3) = uicontrol(handles.MainFigure, 'Style', 'edit',...
                'String', 'Time, s', 'Units', 'characters', ...
                'Position', [controlLeft,controlTop(2),controlWidth,controlHeight], ...
                'Callback', 'kazan(''lDirlist_Callback'',gcbo,[],guidata(gcbo))');
            guidata(handles.MainFigure, handles);
        else % resize
            set(handles.ReadOptHandle(1), 'Position', [controlLeft,controlTop(1),controlWidth2,controlHeight]);
            set(handles.ReadOptHandle(2), 'Position', [controlLeft2,controlTop(1),controlWidth2,controlHeight]);
            set(handles.ReadOptHandle(3), 'Position', [controlLeft,controlTop(2),controlWidth,controlHeight]);
        end
    case  {'*.parElChem_PAR', '*.*ElChem_PAR'} 
        spaceneeded = 2*controlHeight+2*Border;

        % create three controls for unit, delimiter and read format
        if isempty(handles.ReadOptHandle), % create
            handles.ReadOptHandle(1) = uicontrol(handles.MainFigure, 'Style', 'popupmenu',...
                'String', {'I(E)', 'I(T)', 'E(T)'}, ...
                'UserData', {'I(E)', 'I(T)', 'E(T)'}, 'Units', 'characters', ...
                'Position', [controlLeft,controlTop(1),controlWidth,controlHeight], ...
                'Callback', 'kazan(''lDirlist_Callback'',gcbo,[],guidata(gcbo))');
%             handles.ReadOptHandle(2) = uicontrol(handles.MainFigure, 'Style', 'popupmenu',...
%                 'String', {'Y', 'rY, iY', 'X,Y', 'X,rY,iY', 'idx, X, Y', 'idx, X, rY, iY'}, 'Units', 'characters', ...
%                 'Value', 4,'Position', [controlLeft2,controlTop(1),controlWidth2,controlHeight], ...
%                 'Callback', 'kazan(''lDirlist_Callback'',gcbo,[],guidata(gcbo))');
            handles.ReadOptHandle(2) = uicontrol(handles.MainFigure, 'Style', 'edit',...
                'String', '-1', 'Units', 'characters', ...
                'Position', [controlLeft,controlTop(2),controlWidth,controlHeight], ...
                'Callback', 'kazan(''lDirlist_Callback'',gcbo,[],guidata(gcbo))');
            guidata(handles.MainFigure, handles);
        else % resize
            set(handles.ReadOptHandle(1), 'Position', [controlLeft,controlTop(1),controlWidth,controlHeight]);
%             set(handles.ReadOptHandle(2), 'Position', [controlLeft2,controlTop(1),controlWidth2,controlHeight]);
            set(handles.ReadOptHandle(2), 'Position', [controlLeft,controlTop(2),controlWidth,controlHeight]);
        end        
    case {'*.scrHYPERCH'}
        spaceneeded = 2*controlHeight+2*Border;

        % create three controls for unit, delimiter and read format
        if isempty(handles.ReadOptHandle), % create
            handles.ReadOptHandle(1) = uicontrol(handles.MainFigure, 'Style', 'edit',...
                'String', '10', ...
                'Units', 'characters', ...
                'Position', [controlLeft,controlTop(1),controlWidth2,controlHeight], ...
                'Callback', 'kazan(''lDirlist_Callback'',gcbo,[],guidata(gcbo))');
            handles.ReadOptHandle(2) = uicontrol(handles.MainFigure, 'Style', 'edit',...
                'String', '1000', 'Units', 'characters', ...
                'Position', [controlLeft2,controlTop(1),controlWidth2,controlHeight], ...
                'Callback', 'kazan(''lDirlist_Callback'',gcbo,[],guidata(gcbo))');
            handles.ReadOptHandle(3) = uicontrol(handles.MainFigure, 'Style', 'edit',...
                'String', '[1700, 2200]', 'Units', 'characters', ...
                'Position', [controlLeft,controlTop(2),controlWidth,controlHeight], ...
                'Callback', 'kazan(''lDirlist_Callback'',gcbo,[],guidata(gcbo))');
            guidata(handles.MainFigure, handles);
        else % resize
            set(handles.ReadOptHandle(1), 'Position', [controlLeft,controlTop(1),controlWidth2,controlHeight]);
            set(handles.ReadOptHandle(2), 'Position', [controlLeft2,controlTop(1),controlWidth2,controlHeight]);
            set(handles.ReadOptHandle(3), 'Position', [controlLeft,controlTop(2),controlWidth,controlHeight]);
        end        
    case {'*.RMN_2D','*.rmnRMN_2D'}
        spaceneeded = 1*controlHeight+1*Border;
        if isempty(handles.ReadOptHandle), % create
            
            handles.ReadOptHandle(1) = uicontrol(handles.MainFigure, 'Style', 'popupmenu',...
                'String', {'D1:Time D2:Time', 'D1:Time D2:Freq', 'D1:Freq D2:Time', 'D1:Freq D2:Freq'}, ...
                'Units', 'characters', 'UserData', {'2DTT', '2DTF', '2DFT', '2DFF'}, ...
                'Position', [controlLeft,controlTop(1),controlWidth,controlHeight], ...
                'Callback', 'kazan(''lDirlist_Callback'',gcbo,[],guidata(gcbo))');
            guidata(handles.MainFigure, handles);
        else % resize
            set(handles.ReadOptHandle(1), 'Position', [controlLeft,controlTop(1),controlWidth2,controlHeight]);
        end
    case {'*.*AktaStart', '*.resAktaStart'}
        spaceneeded = 1*controlHeight+1*Border;
        if isempty(handles.ReadOptHandle), % create
            
            handles.ReadOptHandle(1) = uicontrol(handles.MainFigure, 'Style', 'popupmenu',...
                'String', {'Volume, mL', 'Time, min'}, ...
                'Units', 'characters', 'UserData', {'Volume, mL', 'Time, min'}, ...
                'Position', [controlLeft,controlTop(1),controlWidth,controlHeight], ...
                'Callback', 'kazan(''lDirlist_Callback'',gcbo,[],guidata(gcbo))');
            guidata(handles.MainFigure, handles);
        else % resize
            set(handles.ReadOptHandle(1), 'Position', [controlLeft,controlTop(1),controlWidth2,controlHeight]);
        end
    case '*.*OPUS_FTIR'
        spaceneeded = 1*controlHeight+1*Border; 
        if isempty(handles.ReadOptHandle), % create
            handles.ReadOptHandle(1) = uicontrol(handles.MainFigure, 'Style', 'popupmenu',...
                'String', {'Interferogram (S_IFG)', 'Data (S_SC)', 'Data-Baseline (_AB)'}, ...
                'Units', 'characters', 'UserData', [1, 2, 3], ...
                'Position', [controlLeft,controlTop(1),controlWidth,controlHeight], ...
                'Callback', 'kazan(''lDirlist_Callback'',gcbo,[],guidata(gcbo))');
            guidata(handles.MainFigure, handles);
        else % resize
            set(handles.ReadOptHandle(1), 'Position', [controlLeft,controlTop(1),controlWidth,controlHeight]);
        end
    case {'*.spaNicoletSPA', '*.*NicoletSPA'}
        spaceneeded = 1*controlHeight+1*Border;
        if isempty(handles.ReadOptHandle), % create
            
            handles.ReadOptHandle(1) = uicontrol(handles.MainFigure, 'Style', 'popupmenu',...
                'String', {'Absorbtion/Transmision', 'Interferogram'}, ...
                'Units', 'characters', 'UserData', {'Absorbtion/Transmision', 'Interferogram'}, ...
                'Position', [controlLeft,controlTop(1),controlWidth,controlHeight], ...
                'Callback', 'kazan(''lDirlist_Callback'',gcbo,[],guidata(gcbo))');
            guidata(handles.MainFigure, handles);
        else % resize
            set(handles.ReadOptHandle(1), 'Position', [controlLeft,controlTop(1),controlWidth2,controlHeight]);
        end        
    case {'*.datASCII2D', '*.*ASCII2D', '*.txtASCII2D', '*.ASCII2D'}
        spaceneeded = 4*controlHeight+4*Border;

        % create three controls for unit, delimiter and read format
        if isempty(handles.ReadOptHandle), % create
            % control for unit
            handles.ReadOptHandle(1) = uicontrol(handles.MainFigure, 'Style', 'edit',...
                'String', 'Time, s', 'Units', 'characters', ...
                'Position', [controlLeft,controlTop(4),sz,controlHeight], ...
                'Callback', 'kazan(''lDirlist_Callback'',gcbo,[],guidata(gcbo))');
            % control for X: label, xsize, dwelltime
            handles.ReadOptHandle(2) = uicontrol(handles.MainFigure, 'Style', 'text',...
                'String', 'X: n,dw', 'Units', 'characters', ...
                'Position', [controlLeft,controlTop(3),controlWidth3,controlHeight],...
                'Callback', 'kazan(''lDirlist_Callback'',gcbo,[],guidata(gcbo))');
            handles.ReadOptHandle(3) = uicontrol(handles.MainFigure, 'Style', 'edit',...
                'String', '1', 'Units', 'characters', ...
                'Position', [controlLeft32,controlTop(3),controlWidth3,controlHeight],...
                'Callback', 'kazan(''lDirlist_Callback'',gcbo,[],guidata(gcbo))');
            handles.ReadOptHandle(4) = uicontrol(handles.MainFigure, 'Style', 'edit',...
                'String', '1', 'Units', 'characters', ...
                'Position', [controlLeft33,controlTop(3),controlWidth3,controlHeight],...
                'Callback', 'kazan(''lDirlist_Callback'',gcbo,[],guidata(gcbo))');
            % control for Y: label, xsize, dwelltime
            handles.ReadOptHandle(5) = uicontrol(handles.MainFigure, 'Style', 'text',...
                'String', 'Y: n,dw', 'Units', 'characters', ...
                'Position', [controlLeft,controlTop(2),controlWidth3,controlHeight],...
                'Callback', 'kazan(''lDirlist_Callback'',gcbo,[],guidata(gcbo))');
            handles.ReadOptHandle(6) = uicontrol(handles.MainFigure, 'Style', 'edit',...
                'String', '1', 'Units', 'characters', ...
                'Position', [controlLeft32,controlTop(2),controlWidth3,controlHeight],...
                'Callback', 'kazan(''lDirlist_Callback'',gcbo,[],guidata(gcbo))');
            handles.ReadOptHandle(7) = uicontrol(handles.MainFigure, 'Style', 'edit',...
                'String', '1', 'Units', 'characters', ...
                'Position', [controlLeft33,controlTop(2),controlWidth3,controlHeight],...
                'Callback', 'kazan(''lDirlist_Callback'',gcbo,[],guidata(gcbo))');
            % control for central frequency
            handles.ReadOptHandle(8) = uicontrol(handles.MainFigure, 'Style', 'text',...
                'String', 'C. Freq', 'Units', 'characters', ...
                'Position', [controlLeft,controlTop(1),controlWidth2,controlHeight],...
                'Callback', 'kazan(''lDirlist_Callback'',gcbo,[],guidata(gcbo))');
            handles.ReadOptHandle(9) = uicontrol(handles.MainFigure, 'Style', 'edit',...
                'String', '0', 'Units', 'characters', ...
                'Position', [controlLeft2,controlTop(1),controlWidth2,controlHeight],...
                'Callback', 'kazan(''lDirlist_Callback'',gcbo,[],guidata(gcbo))');
            guidata(handles.MainFigure, handles);
        else
            set(handles.ReadOptHandle(1), 'Position', [controlLeft,controlTop(4),sz,controlHeight]);
            set(handles.ReadOptHandle(2), 'Position', [controlLeft,controlTop(3),controlWidth3,controlHeight]);
            set(handles.ReadOptHandle(3), 'Position', [controlLeft32,controlTop(3),controlWidth3,controlHeight]);
            set(handles.ReadOptHandle(4), 'Position', [controlLeft33,controlTop(3),controlWidth3,controlHeight]);
            set(handles.ReadOptHandle(5), 'Position', [controlLeft,controlTop(2),controlWidth3,controlHeight]);
            set(handles.ReadOptHandle(6), 'Position', [controlLeft32,controlTop(2),controlWidth3,controlHeight]);
            set(handles.ReadOptHandle(7), 'Position', [controlLeft33,controlTop(2),controlWidth3,controlHeight]);
            set(handles.ReadOptHandle(8), 'Position', [controlLeft,controlTop(1),controlWidth2,controlHeight]);
            set(handles.ReadOptHandle(9), 'Position', [controlLeft2,controlTop(1),controlWidth2,controlHeight]);
        end
    case {'*.*UofC_bin'}
        spaceneeded = 1*controlHeight+1*Border;
        handles.ReadOptHandle(1) = uicontrol(handles.MainFigure, 'Style', 'popupmenu',...
            'String', {'Modula', 'LabView'}, ...
            'Units', 'characters', 'UserData', [1, 2, 3], ...
            'Position', [controlLeft,controlTop(1),controlWidth,controlHeight], ...
            'Callback', 'kazan(''lDirlist_Callback'',gcbo,[],guidata(gcbo))');
        guidata(handles.MainFigure, handles);
    case {'*.*SSRL-EXAFS'}
        spaceneeded = 1*controlHeight+1*Border; 
        if isempty(handles.ReadOptHandle), % create
            handles.ReadOptHandle(1) = uicontrol(handles.MainFigure, 'Style', 'popupmenu',...
                'String', {'SUM(SCA1)/I0', 'I0/RTC', 'LOG(I1/I0)', 'LOG(I1/I2)', 'LOG(I0/I2)',  'SCA1/I0'}, ...
                'Units', 'characters', 'UserData', {'sum(SCA1, 2)./I0', 'I0./rtc', 'log(I1./I0)', 'log(I1./I2)', 'log(I0./I2)',  'SCA1./I0(:, ones(size(SCA1, 2), 1))'}, ...
                'Position', [controlLeft,controlTop(1),controlWidth,controlHeight], ...0
                'Callback', 'kazan(''lDirlist_Callback'',gcbo,[],guidata(gcbo))');
            guidata(handles.MainFigure, handles);
        else % resize
            set(handles.ReadOptHandle(1), 'Position', [controlLeft,controlTop(1),controlWidth,controlHeight]);
        end
    case {'*.*Agilent_BKN', '*.bknAgilent_BKN'}
        spaceneeded = 1*controlHeight+2*Border;
         if isempty(handles.ReadOptHandle), % create
            handles.ReadOptHandle(1) = uicontrol(handles.MainFigure, 'Style', 'popupmenu',...
                'String', {'Dataset #1', 'Dataset #2', 'Dataset #3', 'Dataset #4','Dataset #5','Dataset #6','Dataset #7','Dataset #8','Dataset #9'}, ...
                'Units', 'characters', 'UserData', [1, 2, 3, 4, 5, 6, 7, 8, 9], ...
                'Position', [controlLeft,controlTop(1),controlWidth,controlHeight], ...
                'Callback', 'kazan(''lDirlist_Callback'',gcbo,[],guidata(gcbo))');
            guidata(handles.MainFigure, handles);
        else % resize
            set(handles.ReadOptHandle(1), 'Position', [controlLeft,controlTop(1),controlWidth,controlHeight]);
         end
     case {'*.*ESRxml', '*.xmlESRxml'}
        spaceneeded = 1*controlHeight+2*Border;
         if isempty(handles.ReadOptHandle), % create
            handles.ReadOptHandle(1) = uicontrol(handles.MainFigure, 'Style', 'popupmenu',...
                'String',  {'4096p Reproc', 'Raw', 'Processed'}, ...
                'Units', 'characters', 'UserData', [1, 2], ...
                'Position', [controlLeft,controlTop(1),controlWidth,controlHeight], ...
                'Callback', 'kazan(''lDirlist_Callback'',gcbo,[],guidata(gcbo))');
            guidata(handles.MainFigure, handles);
        else % resize
            set(handles.ReadOptHandle(1), 'Position', [controlLeft,controlTop(1),controlWidth,controlHeight]);
        end
    otherwise
        spaceneeded = 0;
end
% decrease space for directory list
rectD(2) = rectD(2) + spaceneeded;
rectD(4) = rectD(4) - spaceneeded;
set(handles.lDirlist, 'Position', rectD);

MainColor =  get(handles.MainFigure, 'Color');
set(handles.ReadOptHandle, 'BackgroundColor', MainColor);

% % ---------------------------------------------------
% %%%%%%%%%%%%% MATLAB people are idiots, they removed full crosshair !!!!!!!
% function crossHair = local_createCrossHair(fig)
% % Create thin uicontrols with black backgrounds to simulate fullcrosshair pointer.
% % 1: horizontal left, 2: horizontal right, 3: vertical bottom, 4: vertical top
% for k = 1:4
%     crossHair(k) = uicontrol(fig, 'Style', 'text', 'Visible', 'off', 'Units', 'pixels', 'BackgroundColor', [0.5 0.5 0.5], 'HandleVisibility', 'off', 'HitTest', 'off'); %#ok<AGROW>
% end
% end
% function ginputWindowButtonMotion(fig, ~, crossHair)
% % Manages the simulated display of the full cross hair pointer
% updateCrossHair(fig, crossHair);
% end
% function updateCrossHair(fig, crossHair)
% % update cross hair for figure.
% gap = 3; % 3 pixel view port between the crosshairs
% cp = hgconvertunits(fig, [fig.CurrentPoint 0 0], fig.Units, 'pixels', fig);
% cp = cp(1:2);
% figPos = hgconvertunits(fig, fig.Position, fig.Units, 'pixels', fig.Parent);
% figWidth = figPos(3);
% figHeight = figPos(4);
% 
% % Early return if point is outside the figure
% if cp(1) < gap || cp(2) < gap || cp(1)>figWidth-gap || cp(2)>figHeight-gap
%     return
% end
% 
% set(crossHair, 'Visible', 'on');
% thickness = 1; % 1 Pixel thin lines. 
% set(crossHair(1), 'Position', [0 cp(2) cp(1)-gap thickness]);
% set(crossHair(2), 'Position', [cp(1)+gap cp(2) figWidth-cp(1)-gap thickness]);
% set(crossHair(3), 'Position', [cp(1) 0 thickness cp(2)-gap]);
% set(crossHair(4), 'Position', [cp(1) cp(2)+gap thickness figHeight-cp(2)-gap]);
% end


% ---------------------------------------------------
function MainFigure_ButtonDownFcn(h, ev, handles)
switch get(handles.MainFigure, 'SelectionType'),
    case 'open'
        %     if strcmp(get(handles.mToolsZoom, 'Checked'), 'off'), return; end
        handles.Type = 0;
        handles.MMove = 0;
        %     hdl = [handles.mViewAxistight, handles.mViewAxisauto handles.mViewAxisfixed];
        %     hand = findobj(hdl, 'Checked', 'on');
        mViewAxisSet(handles.mViewAxistight, [], handles);
    case 'extend',
        if strcmp(get(handles.mToolsZoom, 'Checked'), 'on'), return; end
        if strcmp(safeget(handles.src.ax, 'type', 'data'), '3dmolecule')   
            handles.Type = 2; % Make orbit camera
            handles.cameraxy = get(handles.MainFigure, 'CurrentPoint');
            handles.cameramax = [diff(get(handles.aReal, 'XLim')), diff(get(handles.aReal, 'YLim'))];
        else
            handles.Type  = 1;
        end
        handles.MMove = 1;
        set(handles.MainFigure, 'Pointer', 'fullcrosshair');
    case 'normal' % for zoom
        if get(handles.MainFigure, 'CurrentObject') == handles.fMoveDataBr % Change width of Directory list
            handles.changePart = 1; % Start change width of Directory list
            set(handles.fMoveDataBr, 'BackgroundColor', [1, 0.9, 0.9]);
        else % zoom
            if strcmp(get(handles.mToolsZoom, 'Checked'), 'off'), return; end
            handles.Type  = 0;
            handles.MMove = 0;
            [x, y, axhand] = getcoord(handles);
            curval = get(handles.MainFigure, 'CurrentPoint');
            if ~isempty(x)
                handles.StartP = [x, y, curval];
                rbbox;
                [x, y] = getcoord(handles, axhand);
                handles.FinishP = [x, y];
                if handles.FinishP(1) == handles.StartP(1) || handles.FinishP(2) == handles.StartP(2), return; end
                if ~isempty(x)
                    set(axhand, 'XLim', sort([handles.StartP(1), handles.FinishP(1)]), ...
                        'YLim', sort([handles.StartP(2), handles.FinishP(2)]));
                    %             if isfield(handles, 'Rectang')
                    %                 delete(handles.Rectang);
                    %
                    %             end
                    % switching zoom off boep maybe good maybe - not
                    % I hope, this is golden middle. alsi
                    if strcmp(handles.autoZoomOff, 'on'),
                        mTools(handles.mToolsZoom,[],handles);
                    end
                    if strcmp(handles.FixedZoomOpt, 'on'),
                        mViewAxisSet(handles.mViewAxisfixed,[],handles);
                    end
                end
            end
        end

end
guidata(handles.MainFigure, handles);
% ---------------------------------------------------
function MainFigure_ButtonUpFcn(h, ev, handles)

if handles.changePart
    MainFigure_ResizeFcn(h, [], handles);
    set(handles.fMoveDataBr, 'BackgroundColor', [0.7, 0.7, 0.7]);
end
handles.changePart = 0;
handles.MMove = 0;

switch handles.Type
    case 1
        set(handles.MainFigure, 'Pointer', 'arrow');
    case 2
        set(handles.MainFigure, 'Pointer', 'arrow');
end
guidata(handles.MainFigure, handles);

% ---------------------------------------------------
function MainFigure_ButtonMotionFcn(h, ev, handles)
% protection
if ~isfield(handles, 'changePart'), return
end
handles.MMove = safeget(handles, 'MMove' , 0);
handles.Type = safeget(handles, 'Type' , 0);

if handles.MMove
    switch handles.Type
        case 1
            [x, y] = getcoord(handles);
            if ~isempty(x)
                s = sprintf(' x= %g, y= %f',x,y);
                set(handles.stCoordinates,'String', s);
            end
        case 2
            aaa = get(handles.MainFigure, 'CurrentPoint');
            phi = (handles.cameraxy(1) - aaa(1))/handles.cameramax(1)*10;
            theta = (handles.cameraxy(2) - aaa(2))/handles.cameramax(2)*10;
            camorbit(phi, theta);
    end
end
if handles.changePart                          % do width of Directory list
    asd = get(handles.MainFigure, 'CurrentPoint');
    mn_pos = get(handles.MainFigure, 'Position');
    MBarPos = get(handles.fMoveDataBr, 'position');
    if ~(asd(2)>mn_pos(4) || asd(2)<0 || asd(1)>mn_pos(3)-30 || asd(1)<15)
        set(handles.fMoveDataBr, 'BackgroundColor', [1, 0.9, 0.9]);
        handles.DataPartSize = asd(1)-MBarPos(3)*0.5; % to make pointer on the middle of the bar
        guidata(handles.MainFigure, handles);
        handles = DataLoadOptions(handles);
    else
        set(handles.fMoveDataBr, 'BackgroundColor', [0.1, 0.1, 0.1]);
    end
end
guidata(handles.MainFigure, handles);
% --------------------------------------------------------------------
function mFavorites(h, eventdata, handles, varargin)
switch get(h, 'Tag')
    case 'mFavorAdd'
        outstr = inimanaging(handles.inifile);
        if ~isempty(outstr),
            if ~isfield(outstr, 'KazanViewer'), outstr.KazanViewer = []; end;
            if isfield(outstr.KazanViewer, 'favorDir')
                outstr.KazanViewer.favorDir{end+1} = get(handles.eDirectory, 'String');
            else
                outstr.KazanViewer.favorDir{1} = get(handles.eDirectory, 'String');
            end
        end
        inimanaging(handles.inifile, outstr);
        CreateFavorites(handles);
    case 'mClearAll'
        menus = get(handles.mFavorites, 'Children');
        delete(findobj(menus, 'UserData', 'favorites'));
        outstr = inimanaging(handles.inifile);
        try outstr.KazanViewer = rmfield(outstr.KazanViewer, 'favorDir'); catch end
        inimanaging(handles.inifile, outstr);
    otherwise
        str = get(h, 'Label');
        set(handles.eDirectory, 'String', str, 'TooltipString', str);
        eDirectory_Callback(handles.eDirectory, [], handles);
end

% --------------------------------------------------------------------
function mFileSave(h, eventdata, handles, varargin)
switch h
    case handles.mFileNewWindow, kazan; return;
    case handles.mFileOpen
        if exist('uigetdir', 'file')
            directory_name = uigetdir(pwd,'Select directory');
        else
            [a, directory_name] = uiputfile([pwd,'\select.dir'],'Select directory');
        end
        if ~directory_name, return; end
        LoadDirectory(directory_name, handles);
    case {handles.mFileSave, handles.mFileSaveProc},
        curdir = pwd;
        cd(handles.dir_path);
        [fname,pname] = uiputfile({'*.dat'; '*.csv';'*.d01';'*.dsc';'*.par'}, 'TYPE THE FILE EXTENSION !!! '); %(don''t forget to type extension !!!)
        cd(curdir);

        if fname == 0, return; end
        [fpath,name,ext] = fileparts(fname);
        src = handles.src; src.fname = handles.fname;

        if h == handles.mFileSaveProc
            [src.y, src.ax] = processing(src.y, src.ax);
            src.ax = getplotaxis(handles.xrepresent, handles.yrepresent, src.ax);
            switch handles.DataPart
                case 1, src.y = real(src.y);
                case 3, src.y = abs(src.y);
            end
        end

        switch lower(ext)
            case '.d01'
                kv_d01write([pname, fname], src);
                disp(['file ',[pname, fname] ,' have been successfully saved in format "d01/exp"'])
            case '.dsc'
                kv_brukerwrite([pname, fname], src);
                disp(['file ',[pname, fname] ,' have been successfully saved in format "dsc/dta"'])
            case '.csv'
                y  = real(src.y);
                y1 = imag(src.y);
                x  = src.ax.x;
                if sum(y1)
                    p = [x,y,y1];
                else
                    p = [x,y];
                end
                csvwrite([pname, fname],p)
                
            case '.par'
                disp('Not implemented')
            otherwise
                if isempty(ext), fname = [fname,'.dat']; end
                kv_asciiwrite([pname, fname], src);
                fid = fopen([pname, fname], 'a');
                fprintf(fid, ';--- created by Kazan Viewer ---\n');
                fprintf(fid,';Original file: %s\n', handles.fname);
                if isfield(src.ax, 'type'), fprintf(fid,';ax.type: %s\n', src.ax.type); end
                if isfield(src.ax, 'xlabel'), fprintf(fid,';ax.xlabel: %s\n', src.ax.xlabel); end
                if isfield(src.ax, 'ylabel'), fprintf(fid,';ax.ylabel: %s\n', src.ax.ylabel); end
%                 if ~isempty(handles.src.ax.filt),
                    if isfield(src.ax, 'filt'), fprintf(fid,';ax.filt: %s\n', handles.src.ax.filt);  end
                    if isfield(src.ax, 'tc'), fprintf(fid,';ax.tc: %f\n', handles.src.ax.tc); end
                    if isfield(src.ax, 'pol'),fprintf(fid,';ax.pol: %d\n', handles.src.ax.pol); end
%                 end
%                 if ~isempty(handles.src.ax.diff),
                    if isfield(src.ax, 'diff'),fprintf(fid,';ax.diff: %s\n', handles.src.ax.diff); end
                    if isfield(src.ax, 'ps_mod'), fprintf(fid,';ax.ps_mod: %f\n', handles.src.ax.ps_mod); end
                    if isfield(src.ax, 'ps_mod_harm'), fprintf(fid,';ax.ps_mod_harm: %f\n', handles.src.ax.ps_mod_harm); end
%                 end
                if isfield(handles.src.ax, 'dx'), fprintf(fid,';ax.dx: %f\n', handles.src.ax.dx); end
                if isfield(handles.src.ax, 'dy'), fprintf(fid,';ax.dy: %f\n', handles.src.ax.dy); end
                if isfield(src.ax, 's'), fprintf(fid,';ax.s: %f\n', src.ax.s); end
                if isfield(handles.src, 'Script') % alsi 28.02.2005
                    fprintf(fid,';Script:\n');
                    for ci = 1:length(handles.src.Script)
                        fprintf(fid,'; %s\n', handles.src.Script{ci});
                    end
                end
                fclose(fid);
                disp(['file ',[pname, fname] ,' successfully saved in format "dat"'])
        end
    case handles.mFileConvert,
        button = questdlg('Continue converting to ASCII (*.dat or *.csv)?', 'Converting all files in the directory to ASCII',...
            'DAT', 'CSV', 'Cancel','DAT');
        if strcmp(button, 'Cancel'), return; end
        for kk=1:length(handles.file_names)
            if handles.is_dir(kk), continue; end
            fname = fullfile(handles.dir_path,handles.file_names{kk});
            hdl = LoadFile(fullfile(handles.dir_path,handles.file_names{kk}), handles);
            y = real(hdl.src.y);
            y1 = imag(hdl.src.y);
            x = hdl.src.ax.x + hdl.src.ax.dx;
            [fpath,name] = fileparts(fname);
           
            if sum(y1)
                p = [x,y,y1];
            else
                p = [x,y];
            end
            if strcmp(button, 'DAT')
                fname = fullfile(fpath, [name,'.dat']);
                save(fname, 'p', '-ascii')
            else
                 fname = fullfile(fpath, [name,'.csv']);
                csvwrite(fname,p)   
            end
            disp([fname,' successfully saved in format "dat"'])
        end
end

% --------------------------------------------------------------------
function mPluginsGlue(h, eventdata, handles, varargin)
posMain = get(handles.MainFigure, 'Position');
kkx = 4;
kky = 9;

try
    for i=1:size(handles.PluginHandle, 2)
        hPlugin = handles.PluginHandle(i);
        posPlugin = get(hPlugin, 'Position');
        switch get(hPlugin, 'Units')
            case 'points',
                posPlugin(2) = (posMain(2) - 2*(i-1) - posPlugin(4)/kky + posMain(4))*kky;
                posPlugin(1) = (posMain(1) + posMain(3) + 1.2 + 8*(i-1))*kkx;
                set(hPlugin, 'Position', posPlugin);
            case 'characters',
                posPlugin(2) = posMain(2) - 2*(i-1) - posPlugin(4) + posMain(4);
                posPlugin(1) = posMain(1) + posMain(3) + 1.2 + 8*(i-1);
                set(hPlugin, 'Position', posPlugin);
        end
    end
catch
    DeletePlugins(hPlugin, handles.MainFigure);
end
guidata(handles.MainFigure, handles);

function mPluginsReload(h, event, handles)
plugs = kazan_plugins;
fpath=handles.pluginsdir;
addpath([fpath, filesep, 'plugins']);

%         plugins = dir(fullfile(handles.pluginsdir, '*.fig'));
callback = 'kazan(''RunPlugins_Callback'',gcbo,[],guidata(gcbo))';

oldMen = get(handles.mPlugins, 'Children');
if ~isempty(oldMen)
    delete(oldMen);
end
%         pnumber = 0;
%         for ii=1:length(plugins)
%             pl = plugins(ii);
%             if strcmp(pl.name, [name ext])==0
%                 [p,n] = fileparts(pl.name);
%                 hmenu = uimenu(handles.mPlugins, 'Label',upper(n), 'Callback', callback, ...
%                     'Tag', ['mPlugins' n], 'UserData', n);
%                 if pnumber < 10
%                     set(hmenu, 'Accelerator', num2str(pnumber));
%                 end
%                 pnumber = pnumber + 1;
%             end
%         end
pnumber = 0;
for ii=1:length(plugs)
    hmenu = uimenu(handles.mPlugins, 'Label', plugs{ii}.name, 'Callback', callback, ...
        'Tag', ['mPlugins' plugs{ii}.name], 'UserData', plugs{ii}.file);
    if pnumber < 10
        set(hmenu, 'Accelerator', num2str(pnumber));
    end
    pnumber = pnumber + 1;
end
handles.mPluginsReposition = uimenu(handles.mPlugins, 'Label','Arrange plugins', 'Tag', 'mPluginsReposition', 'Separator', 'on', 'Callback', ...
    'kazan(''mPluginsGlue'',gcbo,[],guidata(gcbo))');
handles.mPluginsReload = uimenu(handles.mPlugins, 'Label','Reload plugins', 'Tag', 'mPluginsReload', 'Separator', 'of', 'Callback', ...
    'kazan(''mPluginsReload'',gcbo,[],guidata(gcbo))');


% --------------------------------------------------------------------
function RunPlugins_Callback(h, eventdata, handles, varargin)
hPlugIns = eval(get(h, 'UserData'));

% check if plugin exist
if isempty(handles.PluginHandle)
    idx = 0;
else
    idx = handles.PluginHandle == hPlugIns;
end

if sum(idx) == 0
    handles.PluginHandle(end+1)=hPlugIns;
end

gdata = guidata(hPlugIns);
gdata.hh = handles.MainFigure;
guidata(hPlugIns, gdata);
guidata(handles.MainFigure, handles);
[path, name] = fileparts(get(hPlugIns, 'FileName'));
str = [name '(''FileUpdate'',guidata(hPlugIns));'];
try eval(str); catch end


% --------------------------------------------------------------------
function DeletePlugins(phandle, MainFigure)
handles = guidata(MainFigure);
handles.PluginHandle = handles.PluginHandle(handles.PluginHandle~=phandle);
guidata(MainFigure, handles);

% --------------------------------------------------------------------
function res = GetPlotNumber(handles)
if handles.DataSource~=1
    nplots = size(handles.out, 2);
else nplots = 0;
end
% if DataSource 1 or 3 print source
if handles.DataSource~=2,
    fplot = 0;
else fplot = 1;
end
res = fplot:nplots;

% --------------------------------------------------------------------
function PlotData(mainform, varargin)
handles = guidata(mainform);

% this part is used to redirect output to different plots
if nargin < 2, haReal = handles.aReal;
else haReal = varargin{1};
end

if nargin < 3, haImag = handles.aImag;
else haImag = varargin{2};
end

if nargin < 4, fig = handles.MainFigure;
else fig = varargin{3};
end

raxis = axis(haReal);
iaxis = axis(haImag);
newax.xlabel = '?,';
newax.ylabel = '?,';

%  clear all plots
ch = allchild(haReal);
delete(ch);
ch = allchild(haImag);
delete(ch);
set(haReal, 'NextPlot', 'add');
set(haImag, 'NextPlot', 'add');

if strcmp(safeget(handles.src.ax, 'type', 'data'), 'image') && (handles.DataSource==1 || handles.DataSource==3);
    [sxxx, syyy] = size(handles.src.y);
    if floor(syyy/3)==syyy/3
        y = reshape(handles.src.y, [sxxx, syyy/3, 3]);
    end
    image(y, 'Parent', haReal);
    set(haReal, 'YDir', 'reverse')
    axis(haReal, 'image');
elseif strcmp(safeget(handles.src.ax, 'type', 'data'), '3dmolecule')    
   set(handles.MainFigure, 'KeyPressFcn', 'cameratoolbar(''show'');');
% alsi 08.04.2008 plot molecular structure in 3D
    Str = handles.src.ax.cellstruct{1};
    [X, Y, Z]  = sphere(5);
    [cX, cY, cZ] = cylinder(5);
    
    for ii = 1:(Str.natoms-1)
        ratom  = 0.3;
        surf(haReal, X*ratom+Str.XYZ(ii, 1), Y*ratom+Str.XYZ(ii, 2), Z*ratom+Str.XYZ(ii, 3));
        rbond = 0.005;
        bnds = Str.bonds{ii};
        for cc = 1:length(bnds)

            R = Str.XYZ(bnds(cc), :)-Str.XYZ(ii, :);
            d = sqrt(sum(R.^2));

            costh = R(3)/d;
            dprod = sqrt(sum(R(1:2).^2));
            sinth = dprod/d;
            if dprod ==0
                sinph = 0;
                cosph = 1;
            else
                cosph = R(1)/dprod;
                sinph = R(2)/dprod;
            end
            tZ = cZ*d;
            tX = cX*rbond;
            tY = cY*rbond;
            %             nX = costh*cosph*tX + costh*sinph*tY-sinth*tZ;
            %             nY = -sinph*tX + cosph*tY;
            %             nZ = sinth*cosph*tX + sinth*sinph*tY + costh*tZ;

            nX = costh*cosph*tX - sinph*tY + sinth*cosph*tZ;
            nY = costh*sinph*tX + cosph*tY + sinth*sinph*tZ;
            nZ = -sinph*tX + costh*tZ;

            surf(haReal, nX+Str.XYZ(ii, 1), nY+Str.XYZ(ii, 2), nZ+Str.XYZ(ii, 3));
        end
    end
    surf(haReal, X*ratom+Str.XYZ(end, 1), Y*ratom+Str.XYZ(end, 2), Z*ratom+Str.XYZ(end, 3));
    axis equal; 
else
    % protecting axis against upside down etc
    axis(haReal, 'normal');
    axis(haImag, 'normal');
    set(haReal, 'XDir', 'normal')
    set(haReal, 'YDir', 'normal')

    for ii = GetPlotNumber(handles)
        % data selection
        if handles.PlotType~=1, proj = 3;
        else proj = handles.Projection;
        end;
        if ii~=0
            % output
            if isempty(handles.out{ii}.y), continue; end
            tt  = GetSlice(handles.out{ii}, proj, handles.Selection);
            UserData = safeget(handles.out{ii}.ax,'UserData',1);
        else
            % source
            if isempty(handles.src.y), continue; end
            tt  = GetSlice(handles.src, proj, handles.Selection);
            UserData = 1;
        end
        ty = tt.y;
        b = whos('ty');
        clear ty;
        if ~strcmp(b.class, 'uint8')
            [y, ax] = processing(tt.y, tt.ax);
        else
            y = tt.y;
            ax = tt.ax;
        end
        newax = getplotaxis(handles.xrepresent, handles.yrepresent, ax);

        % data printing
        switch handles.PlotType
            case {1,2}
                switch handles.DataPart
                    case 1
                        hhh = superplot(handles.xrepresent, handles.yrepresent, newax, real(y), haReal);
                        set(hhh,'UserData',UserData)
                    case 2
                        hhh = superplot(handles.xrepresent, handles.yrepresent, newax, real(y), haReal);
                        set(hhh,'UserData',UserData)
                        hhh = superplot(handles.xrepresent, handles.yrepresent, newax, imag(y), haImag);
                        set(hhh,'UserData',UserData)
                        axis(haImag, handles.axisset);
                    case 3
                        hhh = superplot(handles.xrepresent, handles.yrepresent, newax, abs(y), haReal);
                        set(hhh,'UserData',UserData)
                end
                if handles.PlotType==2 && isfield(newax,'yslicelabel') && ~isempty(newax.yslicelabel)
                    xmin = min(newax.x) + safeget(newax, 'yslicelabeldx',0);
                    xmax = max(newax.x) + safeget(newax, 'yslicelabeldx',0);
                    ymax = y(end,:) + safeget(newax, 'yslicelabeldy',0);
                    for ii=1:length(newax.y)
                        text((xmax-xmin)*.9+xmin, ymax(ii), sprintf(newax.yslicelabel, newax.y(ii)), 'Parent', haReal);
                    end
                end
            case 3,
                set(0,'CurrentFigure', fig)
                if isfield(ax, 'contour'), cc = ax.contour;
                else cc = .1:.1:.6;
                end

                switch handles.DataPart
                    case 1
                        set(fig, 'CurrentAxes', haReal)
                        maxyy=max(max(real(y))); minyy=min(min(real(y)));
                        if maxyy < 0, sh = minyy; else sh = 0; end
                        contour(newax.x, newax.y, real(y.'), cc*maxyy+sh);
                    case 2
                        set(fig, 'CurrentAxes', haReal);
                        maxyy=max(max(real(y))); minyy=min(min(real(y)));
                        if maxyy < 0, sh = minyy; else sh = 0; end
                        contour(outx, newax.y, real(y.'), cc*(maxyy-minyy)+sh);
                        set(fig, 'CurrentAxes', haImag);
                        maxyy=max(max(imag(y))); minyy=min(min(imag(y)));
                        if maxyy < 0, sh = minyy; else sh = 0; end
                        contour(newax.x, newax.y, imag(y.'), cc*(maxyy-minyy)+sh);
                        axis(haImag, handles.axisset);
                    case 3
                        set(fig, 'CurrentAxes', haReal)
                        maxyy=max(max(abs(y))); 
                        contour(newax.x, newax.y, abs(y.'), cc*maxyy);
                end
            case 4,
                if isfield(ax, 'contour'), cc = ax.contour;
                else cc = .1:.1:.6;
                end
                maxyy=max(max(y));
                minyy=min(min(y));
                if maxyy < 0, sh = minyy;
                else sh = 0;
                end
                
                set(0,'CurrentFigure', fig)
                switch handles.DataPart
                    case 1
                        set(fig, 'CurrentAxes', haReal)
                        clim = cc*real(maxyy)+real(sh);
                        imagesc(newax.x, newax.y, real(y.'), [min(clim), max(clim)]);
                    case 2
                        set(fig, 'CurrentAxes', haReal)
                        clim1 = cc*real(maxyy-minyy)+real(sh);
                        imagesc(newax.x, newax.y, real(y.'), [min(clim1), max(clim1)]);
                        set(fig, 'CurrentAxes', haImag)
                        clim2 = cc*imag(maxyy-minyy)+imag(sh);
                        imagesc(newax.x,newax.y, imag(y.'), [min(clim2), max(clim2)]);
                        axis(haImag, handles.axisset);
                    case 3
                        set(fig, 'CurrentAxes', haReal)
                        clim = cc*(max(max(abs(y)))-min(min(abs(y))));
                        imagesc(newax.x, newax.y, abs(y.'), [min(clim), max(clim)]);
                end
        case 5 % for hyscorebox and image2data <<< AlSi 30.03.05 <<<
            if ii==0|strcmp(safeget(handles.out{ii}, 'forceplot', 'contour'), 'image')  % input data
                if isfield(ax, 'contour'), cc = ax.contour;
                else cc = [.1:.1:.6];
                end
                maxyy=max(max(y));
                minyy=min(min(y));         
                if maxyy < 0, sh = minyy;
                else sh = 0;
                end            
                
                set(0,'CurrentFigure', fig)
                if ~isempty(y)
                    b = whos('y');
                    if strcmp(b.class, 'uint8')
                        
                        if length(size(y))==2
                            [sxxx, syyy] = size(y);
                            if floor(syyy/3)==syyy/3
                                y = reshape(y, [sxxx, syyy/3, 3]);
                            end
                        end
                        image(newax.x, newax.y, y, 'Parent', haReal);
%                         set(haReal, 'YDir', 'reverse')
%                         axis(haReal, 'image');
                    else
                    switch handles.DataPart
                    case 1
                        set(fig, 'CurrentAxes', haReal)
                        clim = cc*real(maxyy-minyy)+real(sh);
                        imagesc(newax.x, newax.y, real(y).', [min(clim), max(clim)]);
                        %               surf(ax.x, ax.y, real(y.'), 'LineStyle', 'none');
                    case 2
                        set(fig, 'CurrentAxes', haReal)
                        clim1 = cc*real(maxyy-minyy)+real(sh);
                        imagesc(newax.x, newax.y, real(y).', [min(clim1), max(clim1)]);
                        %                surf(ax.x, ax.y, real(y.'), 'LineStyle', 'none');
                        set(fig, 'CurrentAxes', haImag)
                        clim2 = cc*imag(maxyy-minyy)+imag(sh);
                        imagesc(newax.x, newax.y, imag(y).', [min(clim2), max(clim2)]);
                        %                surf(ax.x, ax.y, imag(y.'), 'LineStyle', 'none');
                        axis(haImag, handles.axisset);
                    case 3
                        set(fig, 'CurrentAxes', haReal)
                        clim = cc*(max(max(abs(y)))-min(min(abs(y))));
                        imagesc(newax.x, newax.y, abs(y).', [min(clim), max(clim)]);
                        %                surf(ax.x, ax.y, abs(y.'), 'LineStyle', 'none');
                    end
                end
                end
            elseif strcmp(safeget(handles.out{ii}, 'forceplot', 'image'), 'contour')
                set(0,'CurrentFigure', fig);
                if isfield(ax, 'contour'), cc = ax.contour;
                else cc = [.1:.1:.6];
                end
                
                switch handles.DataPart
                    case 1
                        set(fig, 'CurrentAxes', haReal)
                        maxyy=max(max(real(y))); minyy=min(min(real(y)));
                        if maxyy < 0, sh = minyy; else sh = 0; end
                        contour(newax.x, newax.y, real(y.'), cc*maxyy+sh, 'w');
                    case 2
                        set(fig, 'CurrentAxes', haReal);
                        maxyy=max(max(real(y))); minyy=min(min(real(y)));
                        if maxyy < 0, sh = minyy; else sh = 0; end
                        contour(outx, newax.y, real(y.'), cc*(maxyy-minyy)+sh, 'w');
                        set(fig, 'CurrentAxes', haImag);
                        maxyy=max(max(imag(y))); minyy=min(min(imag(y)));
                        if maxyy < 0, sh = minyy; else sh = 0; end
                        contour(newax.x, newax.y, imag(y.'), cc*(maxyy-minyy)+sh, 'w');
                        axis(haImag, handles.axisset);
                    case 3
                        set(fig, 'CurrentAxes', haReal)
                        maxyy=max(max(abs(y))); minyy=min(min(abs(y)));
                        contour(newax.x, newax.y, abs(y.'), cc*maxyy, 'w');
                    end
             else
                switch handles.DataPart
                case 1
                    superplot(handles.xrepresent, handles.yrepresent, ax, ax.y, haReal);
                case 2
                    superplot(handles.xrepresent, handles.yrepresent, ax, ax.y, haReal);
                    superplot(handles.xrepresent, handles.yrepresent, ax, ax.y, haImag);
                    axis(haImag, handles.axisset);
                case 3
                    superplot(handles.xrepresent, handles.yrepresent, ax, ax.y, haReal);
                end
            end
        end 
    end
end
if strcmp(get(handles.mViewAxisfixed, 'Checked'), 'on')
    axis(haReal, raxis);
    axis(haImag, iaxis);
    %             axis(haReal, handles.axisset);
    %             axis(haImag, handles.axisset);
else
    axis(haReal, handles.axisset);
    axis(haImag, handles.axisset);
end
axis(haReal, handles.axisaspect);
axis(haImag, handles.axisaspect);
if strcmp(handles.axisaspect1, 'equal')
    axis(haReal, handles.axisaspect1);
    axis(haImag, handles.axisaspect1);
end
grid(haReal, handles.axisgrid)
grid(haImag, handles.axisgrid)

% set label
set(handles.tUnits, 'String', newax.xlabel);

% --------------------------------------------------------------------
function mViewAxisSet(h, eventdata, handles, varargin)
hdl = [handles.mViewAxistight, handles.mViewAxisauto handles.mViewAxisfixed];
hdl1 = [handles.mViewAxissquare handles.mViewAxisnormal];
handles = guidata(handles.MainFigure);
switch h
    case handles.mViewAxistight
        axis(handles.aReal, 'tight')
        axis(handles.aImag, 'tight')
        handles.axisset = 'tight';
    case handles.mViewAxisauto
        axis(handles.aReal, 'auto')
        axis(handles.aImag, 'auto')
        handles.axisset = 'auto';
    case handles.mViewAxisfixed
        handles.axisset = 'auto';
    case handles.mViewAxisnormal
        axis(handles.aReal, 'normal')
        axis(handles.aImag, 'normal')
        handles.axisaspect = 'normal';
        set(h,'Checked','on');
        set(hdl1(hdl1~=h),'Checked','off');
        guidata(handles.MainFigure, handles);
        return
    case handles.mViewAxissquare
        axis(handles.aReal, 'square')
        axis(handles.aImag, 'square')
        handles.axisaspect = 'square';
        set(h,'Checked','on');
        set(hdl1(hdl1~=h),'Checked','off');
        guidata(handles.MainFigure, handles);
        return
    case handles.mViewAxisEqual
        if strcmp(get(h,'Checked'),'on')
            handles.axisaspect1 = 'normal';
            set(h,'Checked','off')
        else
            handles.axisaspect1 = 'equal';
            set(h,'Checked','on')
        end
        axis(handles.aReal, handles.axisaspect1)
        axis(handles.aImag, handles.axisaspect1)
        guidata(handles.MainFigure, handles);
        return
    case handles.mViewShowGrid
        if strcmp(handles.axisgrid,'on'), handles.axisgrid = 'off';
        else handles.axisgrid = 'on';
        end
        grid(handles.aReal, handles.axisgrid)
        grid(handles.aImag, handles.axisgrid)
        guidata(handles.MainFigure, handles);
        set(h,'Checked',handles.axisgrid);
        return
end
set(h,'Checked','on');
set(hdl(hdl~=h),'Checked','off');
guidata(handles.MainFigure, handles);

% --------------------------------------------------------------------
function DataSource_Callback(h, eventdata, handles, varargin)

switch h
    case {handles.rbSource, handles.mViewSource}
        SetDataSource(1, handles)
    case {handles.rbOutput, handles.mViewOutput}
        SetDataSource(2, handles)
    otherwise
        SetDataSource(3, handles)
end

% --------------------------------------------------------------------
function SetDataSource(source, handles)
switch source
    case 1
        set(handles.rbOutput,'Value',0);
        set(handles.rbSource,'Value',1);
        set(handles.mViewOutput,'Checked','off');
        set(handles.mViewSource,'Checked','on');
        set(handles.mViewSnO,'Checked','off');
        handles.DataSource = 1;
    case 2
        set(handles.rbOutput,'Value',1);
        set(handles.rbSource,'Value',0);
        set(handles.mViewOutput,'Checked','on');
        set(handles.mViewSource,'Checked','off');
        set(handles.mViewSnO,'Checked','off');
        handles.DataSource = 2;
    case 3
        set(handles.rbOutput,'Value',1);
        set(handles.rbSource,'Value',1);
        set(handles.mViewOutput,'Checked','on');
        set(handles.mViewSource,'Checked','on');
        set(handles.mViewSnO,'Checked','on');
        handles.DataSource = 3;
    otherwise
        set(handles.rbOutput,'Value',0);
        set(handles.rbSource,'Value',0);
        set(handles.mViewOutput,'Checked','off');
        set(handles.mViewSource,'Checked','off');
        handles.DataSource = 3;
end
guidata(handles.MainFigure, handles);
UpdateSlider(handles);

% --------------------------------------------------------------------
function DataPart_Callback(h, eventdata, handles, varargin)
switch h
    case {handles.rbReal, handles.mViewReal}
        SetDataPart(1, handles);
    case {handles.rbImag, handles.mViewImag}
        SetDataPart(2, handles);
    otherwise
        SetDataPart(3, handles);
end

% --------------------------------------------------------------------
function PlotType_Callback(h, eventdata, handles, varargin)
switch h
    case {handles.tb1D, handles.mView1D}
        SetPlotType(1, handles);
    case {handles.tbStacked, handles.mViewStacked}
        SetPlotType(2, handles);
    case {handles.tbContour, handles.mViewContour}
        SetPlotType(3, handles);
    otherwise
        SetPlotType(4, handles);
end
% --------------------------------------------------------------------
function SetDataPart(type, handles)
hdl = [handles.rbReal, handles.rbImag, handles.rbMagnitude];
hdl1 = [handles.mViewReal, handles.mViewImag, handles.mViewMagnitude];
switch type
    case 1
        set(handles.rbReal,'Value',1);
        set(hdl(hdl ~= handles.rbReal),'Value',0);
        set(handles.mViewReal,'Checked','on');
        set(hdl1(hdl1~=handles.mViewReal),'Checked','off');
        handles.DataPart = 1;
    case 2
        set(handles.rbImag,'Value',1);
        set(hdl(hdl ~= handles.rbImag),'Value',0);
        set(handles.mViewImag,'Checked','on');
        set(hdl1(hdl1~=handles.mViewImag),'Checked','off');
        handles.DataPart = 2;
    otherwise
        set(handles.rbMagnitude,'Value',1);
        set(hdl(hdl ~= handles.rbMagnitude),'Value',0);
        set(handles.mViewMagnitude,'Checked','on');
        set(hdl1(hdl1~=handles.mViewMagnitude),'Checked','off');
        handles.DataPart = 3;
end

guidata(handles.MainFigure, handles);
MainFigure_ResizeFcn(0, 0, handles);
PlotData(handles.MainFigure);

% --------------------------------------------------------------------
function SetPlotType(type, handles)
hdl = [handles.tb1D,  handles.tbStacked, handles.tbContour, handles.tbDensity];
hdl1 = [handles.mView1D, handles.mViewStacked, handles.mViewContour, handles.mViewDensity];

isScroller = handles.isScroller;
handles.isScroller = (type == 1);

switch type
    case 1
        set(handles.tb1D,'Value',1);
        set(hdl(hdl ~= handles.tb1D),'Value',0);
        set(handles.mView1D,'Checked','on');
        set(hdl1(hdl1~=handles.mView1D),'Checked','off');
        handles.PlotType = 1;
    case 2
        set(handles.tbStacked,'Value',1);
        set(hdl(hdl ~= handles.tbStacked),'Value',0);
        set(handles.mViewStacked,'Checked','on');
        set(hdl1(hdl1~=handles.mViewStacked),'Checked','off');
        handles.PlotType = 2;
    case 3
        set(handles.tbContour,'Value',1);
        set(hdl(hdl ~= handles.tbContour),'Value',0);
        set(handles.mViewContour,'Checked','on');
        set(hdl1(hdl1~=handles.mViewContour),'Checked','off');
        handles.PlotType = 3;
    otherwise
        set(handles.tbDensity,'Value',1);
        set(hdl(hdl ~= handles.tbDensity),'Value',0);
        set(handles.mViewDensity,'Checked','on');
        set(hdl1(hdl1~=handles.mViewDensity),'Checked','off');
        handles.PlotType = 4;
end

if isScroller ~= handles.isScroller
    str = {'off', 'on'};
    MainFigure_ResizeFcn([], [], handles);
    set(handles.pb2DSelect, 'Visible', str{handles.isScroller + 1});
    set(handles.sl2Dplot, 'Visible', str{handles.isScroller + 1});
end

% incorrect ...
if handles.isScroller,
    dim = size(handles.src.y);
    if dim(2) > 1
        set(handles.pb2DSelect, 'String', 'X');
        set(handles.sl2Dplot, 'Min', 1, 'Max', dim(2), 'Value', 1,...
            'SliderStep', [1 10]*1/dim(2));
        handles = sl2Dplot_Callback(handles.sl2Dplot, [], handles);
    end
else
    set(handles.stCoordinates, 'String', '');
end

guidata(handles.MainFigure, handles);
PlotData(handles.MainFigure);
% --------------------------------------------------------------------
function handles = pmDatatype_Callback(h, eventdata, handles, varargin)
sel = get(handles.pmDatatype, 'Value');
if sel < 1
    return
end
str = get(handles.pmDatatype, 'String');
handles.Datatype = str{sel};
handles = DataLoadOptions(handles);
guidata(handles.MainFigure, handles);

% --------------------------------------------------------------------
function handles = LoadDirectory(dir_path, handles)
last_pwd = pwd;
if ispc
    if size(dir_path,2) > 0
        if sum(dir_path == ':')==0
            try cd (fullfile(handles.dir_path, dir_path)); catch end
        else
            cd (dir_path);
        end
    else
        try cd (handles.dir_path); catch end
    end
else
    try 
        cd (fullfile(handles.dir_path, dir_path)); 
    catch
        try
            cd (dir_path);
        catch
            disp(['Cannot get to : "', dir_path,'"']);    
        end
    end
end
handles.dir_path = pwd;
dir_dir_struct = dir(handles.dir_path);
dir_dir_struct = dir_dir_struct([dir_dir_struct.isdir]);
if sum(findstr([dir_dir_struct.name], '.'))==6
    dir_dir_struct = dir_dir_struct(3:end);
end
[sorted_dirnames,sorted_dirindex] = sortrows({dir_dir_struct.name}');

dir_size = size(dir_dir_struct, 1);

if ispc
    dir_struct = dir(fullfile(handles.dir_path, handles.ext));
else
    dir_struct = dir(fullfile(handles.dir_path, lower(handles.ext)));
    dir_struct_tmp = dir(fullfile(handles.dir_path, upper(handles.ext)));
    for ii=1:length(dir_struct_tmp), 
        dir_struct(end+1,1)=dir_struct_tmp(ii); 
    end
end
dir_struct = dir_struct(~[dir_struct.isdir]);
file_size = size(dir_struct, 1);

[sorted_names,sorted_index] = sortrows({dir_struct.name}');
sorted_index = sorted_index+dir_size;

% final_names =
final_index = [sorted_dirindex; sorted_index];
final_names = {};
% str{1} = '<HTML><FONT color="blue">';
% str{2} = '<HTML><FONT color="red">';
% str{3} = '<HTML><FONT color="green">';
for ii=1:dir_size
%     final_names{ii} = ['<HTML><FONT color="blue"> [', sorted_dirnames{ii}, ']</FONT></HTML>'];
    final_names{ii} = [handles.DirStartCode, sorted_dirnames{ii},  handles.DirEndCode];
%     [handles.DirStartCode ,aa, handles.DirEndCode]
end
for ii=1:file_size
    final_names{ii+dir_size} = sorted_names{ii};
end

handles.file_names = final_names;
handles.is_dir = [ones(1, dir_size) dir_struct.isdir];
handles.sorted_index = [final_index];
guidata(handles.MainFigure,handles);
set(handles.lDirlist,'String',handles.file_names,'Value',1)
set(handles.eDirectory,'String',pwd, 'TooltipString', pwd)
cd(last_pwd);

% --------------------------------------------------------------------
function lDirlist_Callback(h, eventdata, handles, varargin)

index_selected = get(handles.lDirlist,'Value');
file_list = get(handles.lDirlist,'String');
if isempty(index_selected) return; end
filename = file_list{index_selected}; % Item selected in list box
if handles.is_dir(handles.sorted_index(index_selected)) % If directory
    %strip directory brackets
    nst = numel(handles.DirStartCode)+1;
    nen = numel(filename)-numel(handles.DirEndCode);
    filename = filename(nst:nen);
    
    % Load list box with new directory when double clicked
    seltype = get(handles.MainFigure, 'SelectionType');
    if strcmp(seltype, 'open'), handles = LoadDirectory(filename, handles); end
else
    handles.src.ax = [];
    handles.src.y = [];
    
    handles=LoadFile(fullfile(handles.dir_path, filename), handles);

    smSelPar_Callback([], [], handles);
    guidata(handles.MainFigure,handles);

    % call update and handles repair for all plugins
    if ~isempty(handles.PluginHandle)
        for k=1:size(handles.PluginHandle, 2)
            try
                [pathstr,name] = fileparts(get(handles.PluginHandle(k), 'FileName'));
                str = [name, '(''FileUpdate'',guidata(handles.PluginHandle(k)));'];
                try eval(str); catch end
            catch DeletePlugins(handles.PluginHandle(k), handles.MainFigure);
            end
        end
    end
    handles = guidata(handles.MainFigure);

    source = handles.DataSource;
    guidata(handles.MainFigure, handles);
    if source ~=1, source = 3; end;
    SetDataSource(source, handles);
end

% --------------------------------------------------------------------
function eDirectory_Callback(h, eventdata, handles, varargin)
d = get(handles.eDirectory, 'String');
if size(d, 2)==1
    if ispc
        d = [d ':\'];
    end
end
LoadDirectory(d, handles);

% --------------------------------------------------------------------
function pbUp_Callback(h, eventdata, handles, varargin)
d= handles.dir_path;
str = d;
while ~isempty(str)
    if ispc
        [aa, str] = strtok(str, '\');
    else
        [aa, str] = strtok(str, '/');
    end
end
LoadDirectory('..', handles);
if ~isempty(aa)
    dirs = get(handles.lDirlist, 'String');
    whichone = strcmp(dirs, [handles.DirStartCode ,aa, handles.DirEndCode]);
    k = find(whichone);
    set(handles.lDirlist, 'Value', k);
end

% --------------------------------------------------------------------
function handles = LoadFile(filename, handles)
fileinfo = dir(filename);


handles.lscript = {};
if get(handles.tb1D, 'Value')
    cmContexMenu2(handles.tb1D, [], handles);
end
datatype = [handles.ext handles.Datatype];
handles.fname = filename;
show_mess = 0;
try
    switch datatype
        case {'*.*XWINNMR', '*.XWINNMR'}
            handles.lscript{1} = ['[ax,y,dsc]=XWINNMR(''', filename, ''');'];
        case {'*.*SSRL-EXAFS'}
            num = get(handles.ReadOptHandle(1), 'Value');
            str = get(handles.ReadOptHandle(1), 'UserData');
            handles.lscript{1} = ['[ax,y,dsc]=ssrlexafs(''', filename, ''', ''',str{num} , ''');'];
        case {'*.*WINNMR', '*.WINNMR', '*.aqsWINNMR'}
            handles.lscript{1} = ['[ax,y,dsc]=WINNMR(''', filename, ''');'];
        case {'*.*WEBMOSS'}
            handles.lscript{1} = ['[ax,y,dsc]=moss_webbin(''', filename, ''');'];
        case {'*.rmnRMN_1D', '*.*RMN_1D', '*.RMN_1D'}
            handles.lscript{1} = ['[ax,y,dsc]=rmnread(''', filename, ''', ''type'', ''1D'');'];
        case {'*.rmnRMN_2D', '*.RMN_2D','*.*RMN_2D'}
            idx = get(handles.ReadOptHandle(1), 'Value');
            type = get(handles.ReadOptHandle(1), 'UserData');
            type = type{idx};
            handles.lscript{1} = ['[ax,y,dsc]=rmnread(''', filename, ''', ''type'', ''2D'',''subtype'',''',type,''');'];
        case {'*.datRMNdat'}
            show_mess = 1; 
            handles.lscript{1} = ['[ax,y,dsc]=rmnread(''', filename, ''', ''type'', ''1dtxt'');'];
        case {'*.datASCII2D' '*.*ASCII2D', '*.txtASCII2D', '*.ASCII2D'}
            xsize = get(handles.ReadOptHandle(3), 'String');
            ysize = get(handles.ReadOptHandle(6), 'String');
            xdw = get(handles.ReadOptHandle(4), 'String');
            ydw = get(handles.ReadOptHandle(7), 'String');
            cfreq = get(handles.ReadOptHandle(9), 'String');
            handles.lscript{1} = ['[ax,y,dsc]=rmnread(''', filename, ''', ''type'', ''2dtxt'',', ...
                '''xsize'',',xsize,',''xdwell'',',xdw,',''ysize'',',ysize,',''ydwell'',',ydw,',''cfreq'',',cfreq,');'];
            handles.lscript{end+1} = ['ax.xlabel=''',get(handles.ReadOptHandle(1), 'String'),''';'];
            handles.lscript{end+1} = ['ax.ylabel=''',get(handles.ReadOptHandle(1), 'String'),''';'];
        case {'*.dscXEPR', '*.parXEPR'}
            handles.lscript{1} = ['[ax,y,dsc]=brukerread(''', filename, ''');'];
        case {'*.dscXEPR_JSS', '*.parXEPR_JSS'}
            handles.lscript{1} = ['[ax,y,dsc]=brukerreadJSS(''', filename, ''');'];
        case  {'*.expWIS', '*.d00WIS'}
            handles.lscript{1} = ['[ax,y,dsc]=d00read(''', filename, ''');'];
        case  {'*.dattransient'}
            handles.lscript{1} = ['[ax,y,dsc]=ESPtransread(''', filename, ''');'];
        case  {'*.expSpecMan', '*.d01SpecMan'}
            handles.lscript{1} = ['[ax,y,dsc]=kv_d01read(''', filename, ''');'];
        case  '*.*OPUS_FTIR'
            nn = get(handles.ReadOptHandle(1), 'Value');
            handles.lscript{1} = ['[ax,y,dsc]=opus_read(''', filename, ''', ', num2str(nn) ,');'];
        case  '*.cmblVernierCMBL'
            handles.lscript{1} = ['[ax,y,dsc]=kv_cmblread(''', filename, ''');'];
        case  {'*.spPerkinIR', '*.*PerkinIR'}
            handles.lscript{1} = ['[ax,y,dsc]=perkinelmer_read(''', filename, ''');'];       
        case  {'*.spaNicoletSPA', '*.*NicoletSPA'}
            num = get(handles.ReadOptHandle(1), 'Value');
            handles.lscript{1} = ['[ax,y,dsc]=SPAread(''', filename, ''', ', num2str(num) ,');'];       
        case  '*.scrHYPERCH'
            lw = get(handles.ReadOptHandle(1), 'String');
            npoints = get(handles.ReadOptHandle(2), 'String');
            range = get(handles.ReadOptHandle(3), 'String');
            handles.lscript{1} = ['[ax,y,dsc]=hyperchemread(''', filename, ''', ', lw ,', ', npoints, ', ', range ');'];            
        case {'*.datASCII', '*.*ASCII', '*.ASCII', '*.txtASCII'}
            show_mess = 1;
            num = get(handles.ReadOptHandle(1), 'Value');
            str = get(handles.ReadOptHandle(1), 'UserData');
            delimiter = [str{num}];
            fstcol = get(handles.ReadOptHandle(2), 'Value');
            dunit = get(handles.ReadOptHandle(3), 'String');
            handles.lscript{1} = ['[ax,y,dsc]=asciiread(''', filename, ''',''',delimiter,''',',num2str(fstcol>2),');'];
            if fstcol==5 || fstcol==6
                handles.lscript{end+1} = 'ax.x = y(:,1); y = y(:,2:end);';
            end
            switch fstcol
                case {2,4,5,6},
                    handles.lscript{end+1} = 'ysize = size(y,2);';
                    handles.lscript{end+1} = 'if ~mod(ysize, 2), yy1 = reshape(y, size(y, 1), size(y,2)/2, 2); end;';
                    handles.lscript{end+1} = 'if ~mod(ysize, 2), y = yy1(:,:,1)+i*yy1(:,:,2); end;';
            end
            handles.lscript{end+1} = ['if strcmp(ax.xlabel, ''?''), ax.xlabel=''',dunit,''';end;'];
        case '*.TecMag'
            handles.lscript{1} = ['[ax,y,dsc]=TecMagread(''', filename, ''');'];
        case '*.Spinsight'
            handles.lscript{1} = ['[ax,y,dsc]=kv_read_spinsight(''', filename, ''');'];
        case {'*.jpgimage', '*.*image'}
            handles.lscript{1} = ['y = imread(''', filename, ''');'];
            handles.lscript{end+1} = 'ax.type=''image'';';
            handles.lscript{end+1} = 'ax.xlabel = ''Image width, pixels'';';
            handles.lscript{end+1} = 'dsc = struct(''comment'', ''no info'');';
        case {'*.*UofC_bin'}
            num = get(handles.ReadOptHandle(1), 'Value');
            handles.lscript{1} = ['[ax,y,dsc] = kv_read_halpern(''', filename, ''',''format'',',num2str(num),');'];
        case {'*.*AktaStart', '*.resAktaStart'}
            str = get(handles.ReadOptHandle(1), 'String');
            val = get(handles.ReadOptHandle(1), 'Value');
            handles.lscript{1} = ['[ax,y,dsc] = kv_read_RESAktaStart(''', filename, ''',''',str{val},''');'];
        case {'*.*SXYZ', '*.xyzSXYZ'}
            handles.lscript{1} = ['[ax,dsc] = sxyzread(''', filename, '''); y = 1;'];
        case {'*.*Wband_Brln', '*.wcwWband_Brln', '*.w2pWband_Brln', '*.wtrWband_Brln'}
            handles.lscript{1} = ['[ax,y,dsc]=kv_readWband(''', filename, ''');'];          
        case {'*.parElChem_PAR', '*.*ElChem_PAR'} 
            num = get(handles.ReadOptHandle(1), 'Value');
            str = get(handles.ReadOptHandle(1), 'UserData');
%             segment = eval(str);
% %             fstcol = get(handles.ReadOptHandle(2), 'Value');
            segment = get(handles.ReadOptHandle(2), 'String');
            handles.lscript{1} = ['[ax,y,dsc]=parread(''', filename, ''',',segment,',''',str{num},''');'];
        case {'*.mosMOSS', '*.*MOSS'}
            handles.lscript{1} = ['[ax,y,dsc]=kv_mossread(''', filename, ''');'];
        case {'*.datCST', '*.txtCST', '*.*CST'}    
            handles.lscript{1} = ['[ax,y,dsc]=kv_cstdata(''', filename, ''');'];
        case {'*.*Agilent_BKN', '*.bknAgilent_BKN'}
            nn = get(handles.ReadOptHandle(1), 'Value');
            handles.lscript{1} = ['[ax,y,dsc]=bkn_read(''', filename, ''', ', num2str(nn) ,');'];   
        case {'*.*ESRxml', '*.xmlESRxml'}
            nn = get(handles.ReadOptHandle(1), 'Value');
            rr = 'srp'; % p=processed, r=raw, s=reprocessed
            handles.lscript{1} = ['[ax,y,dsc]=ESRxmlread(''', filename, ''', ''', rr(nn) ,''');'];               
        otherwise
            error('Unknown format');
    end
    % load script:
    for k = 1:size(handles.lscript, 2)
        eval(handles.lscript{k});
    end
    handles.lscript{end+1} = '% plot(ax.x,real(y),ax.x,imag(y))';
    
    if (fileinfo.bytes >= 2^20)&show_mess,
    button = questdlg('Do you want to continue reading?',...
        ['File size ' num2str(fileinfo.bytes/2^20, '%10.2f'), 'Mb'],...
        'Yes','No','Yes');
    if strcmp(button, 'No'), return; end
end
catch
    set(handles.MainFigure, 'Name', ['KAZAN - Data load error']);
    error(['''',handles.Datatype,''' Data load error (probably wrong type). ']);
    ax = struct('x',[], 'xlabel', '?', 'ylabel', '?');
    y  = [];  dsc = [];
end

% this will convert data into 2 dimensional variant
sz = size(y);
if size(sz,2) > 2
  warning('Data dimension is reduced.');  
  y = reshape(y, [sz(1), prod(sz(2:end))]);
  ax.y = [1:prod(sz(2:end))].';
  sz = size(y);
end

% first 'ax' dimension, protection
if ~isfield(ax, 'x')
    ax.x = (1:sz(1)).';
end
if ~isfield(ax, 'xlabel')
    ax.xlabel = 'unknown, pnt';
end
% second 'ax' dimension, protection
if ~isfield(ax, 'y')
    ax.y = [1:size(y, 2)].';
end
if ~isfield(ax, 'ylabel')
    ax.ylabel = 'unknown, pnt';
end
handles.src.y   = y;
handles.src.ax  = ax;
handles.src.dsc = dsc;
handles.src.filename = filename;

handles.src.ax.ext = handles.ext;
handles.src.ax.Datatype = handles.Datatype;
handles.src.ax.filt = handles.filter;
handles.src.ax.diff = handles.process;

% load parameters
str = get(handles.pmSelPar, 'String');
for k=1:size(str, 1)
    fname = str{k};
    if ~isfield(ax, fname)
        if handles.prop.types{k}=='s'
            if ~isempty(handles.prop.values{k})
                eval(['tmp =''', handles.prop.values{k}(handles.prop.values{k}~=''''), ''';']);
            else tmp = '';
            end
        else
            eval(['tmp =', handles.prop.values{k}, ';']);
        end
        handles.src.ax = setfield(handles.src.ax, fname, tmp);
    else
        handles.prop.values{k} = num2str(getfield(handles.src.ax, fname));
    end
end

UpdateSlider(handles);
[PATHSTR,NAME,EXT] = fileparts(filename);

if isfield(handles.src.ax, 'title')
    set(handles.MainFigure, 'Name', ['KAZAN - ',handles.src.ax.title, ' [',NAME,']']);
else
    set(handles.MainFigure, 'Name', ['KAZAN - ',handles.process, ' [',NAME,']']);
end

% --------------------------------------------------------
function mTools(h, eventdata, handles, varargin)
switch h
    case handles.mToolsClearSource
        handles.src.y=[];
        handles.src.ax.x=[];
        handles.src.ax.y=[];
        handles.src.dsc=[];
        guidata(handles.MainFigure, handles);
        PlotData(handles.MainFigure);
        return;
    case handles.mToolsFixProcessing
        handles.src.ax.filter   = safeget(handles,'filter','');
        handles.src.ax.filtpar1 = safeget(handles,'filterpar1','auto');
        handles.src.y=processing(handles.src.y,handles.src.ax);
        handles=cmContexMenu1(handles.cmNormal,[],handles);
        handles=cmContexMenu2(handles.cmxy,[],handles);
        handles=Filter_CallBack(handles.cmNofilter,[],handles);
        guidata(handles.MainFigure, handles);
        PlotData(handles.MainFigure);
        return;
    case handles.mToolsFixSlice
        if handles.PlotType == 1
            handles.src = GetSlice(handles.src,handles.Projection,handles.Selection);
            handles.Projection = 1;
            handles.AntiProjection = 2;
            guidata(handles.MainFigure, handles);
            UpdateSlider(handles);
        end
        return;
    case handles.mToolsLoadscript
        f = '';
        dim = size(handles.lscript, 2);
        for i=1:dim; f = [f sprintf('\n') handles.lscript{1,i}]; end
        disp('Data load script is copied to the clipboard.');
        clipboard('copy',f)
        disp(f);
        return;
    case handles.mToolsCopyLoadscript
        f = '';
        dim = size(handles.lscript, 2);
        for i=1:dim; f = [f sprintf('\n') handles.lscript{1,i}]; end
        disp('Data load script is copied to the clipboard.');
        clipboard('copy',f)
%         disp(f);
        return;
    case handles.mToolsExecLoadscript
        dim = size(handles.lscript, 2);
        for i=1:dim,
            disp(handles.lscript{1,i});
            evalin('base', handles.lscript{1,i}, '');
        end
        return;
    case handles.mToolsZoom
        state = get(handles.mToolsZoom, 'Checked');
        if strcmp('off', state),
            val = 'on';
            P(1:16, 1:16) = NaN;
            P([1, 2, 4:16], 3)= 1;
            P(3, [1, 2, 4:16])= 1;
            P([4, 14], 7:11) = 1;
            P(7:11, [4, 14]) = 1;
            P([6, 12], [5, 13]) = 1;
            P([5, 13], [6, 12]) = 1;
            P(13, 13) = 1;
            P([13, 14], [14, 13]) = 1;
            P([15, 16], [16, 15]) = 1;
            P(7, 11) = 1;
            P([8, 9, 10], 12) = 1;
            P(11, 11) = 1;
            set(handles.MainFigure,'Pointer','custom','PointerShapeCData',P,...
                'PointerShapeHotSpot', [3, 3]);
        else
            val = 'off';
            set(handles.MainFigure,'Pointer','arrow');
        end
        set(handles.mToolsZoom, 'Checked', val);
        return;
    case handles.mToolsAutoZoomOff
        str = {'off', 'on'};
        stat = get(handles.mToolsAutoZoomOff, 'Checked');
        idx = find(strcmp(str, stat));
        set(handles.mToolsAutoZoomOff, 'Checked', str{~(idx-1) + 1});
        handles.autoZoomOff = str{~(idx-1) + 1};
        guidata(handles.MainFigure, handles);
        return;
    case handles.mToolsZoomFixed
        str = {'off', 'on'};
        stat = get(handles.mToolsZoomFixed, 'Checked');
        idx = find(strcmp(str, stat));
        set(handles.mToolsZoomFixed, 'Checked', str{~(idx-1) + 1});
        handles.FixedZoomOpt = str{~(idx-1) + 1};
        guidata(handles.MainFigure, handles);
        return;
    case {handles.mCoord, handles.tbMeasure}
        p=[];
        disp('press any key to cancel')
        while(1)
            [x,y, button]=ginput(1);
            if isempty(button), button=2; end;
            switch button
                case 1,
                    if handles.DataSource==2
                        newax = getplotaxis(handles.xrepresent, handles.yrepresent, handles.out{1}.ax);
                    else
                        newax = getplotaxis(handles.xrepresent, handles.yrepresent, handles.src.ax);
                    end
                    if handles.Projection == 1 || handles.PlotType > 2, xx = newax.x; yy = newax.y; else xx = newax.y; yy = newax.x; end

                    [vvv,xpnt] = min(abs(xx-x));
                    if handles.PlotType < 2
                        s = sprintf('[%d] x= %g, y= %g',xpnt, x,y);
                    elseif handles.PlotType == 2
                        ddyy = handles.src.y(xpnt, :);
                        [vvv,ypnt] = min(abs(ddyy-y));
                        s = sprintf('[%d,%d] x= %f, y= %f', xpnt, ypnt, x,y);
                    else
                        [vvv,ypnt] = min(abs(yy-y));
                        s = sprintf('[%d,%d] x= %f, y= %f', xpnt, ypnt, x,y);
                    end
                    set(handles.stCoordinates,'String', s);
                    disp(s)
                    p(end+1, :) = [x,y];
                otherwise
                    break;
            end
        end;
        switch size(p,1)
            case 2,
                % distance
                s = sprintf('dist x= %g,    y= %f',p(2,1)-p(1,1),p(2,2)-p(1,2));
                disp(s)
                set(handles.stCoordinates,'String', s);
                disp(sprintf('center x= %g,    y= %f',(p(2,1)+p(1,1))/2,(p(2,2)+p(1,2))/2));
            case 4,
                % amplitude relation
                s = sprintf('(y4-y3)/(y2-y1) = %f',(p(4,2)-p(3,2))/(p(2,2)-p(1,2)));
                set(handles.stCoordinates,'String', s);
                disp(sprintf('(x4-x3)/(x2-x1) = %f',(p(4,1)-p(3,1))/(p(2,1)-p(1,1))));
                disp(s)
        end
        disp(sprintf('\n'))
    case handles.mToolsAssociatedWindowReset
        AssociateFigure(handles.MainFigure,[]);
    case handles.mToolsShowinwindow0
        a = figure;
        set(a,'UserData', struct('MainFigure',handles.MainFigure));
        pMenu = uimenu('Parent', a, 'Label', 'KV_Tools');
        uimenu('Parent', pMenu, 'Label', 'Export Image', 'Callback', 'kv_exportimage(gcf, ''gui'')');  % alsi 12.01.05
        uimenu('Parent', pMenu, 'Label', 'Profile Manager', 'Callback', 'kv_printprofile');  % alsi 12.01.05
        uimenu('Parent', pMenu, 'Label', 'Associate',...
            'Callback', ['ud = get(gcf ,''UserData''); hand = getfield(ud,''MainFigure''); gdata = guidata(hand);',...
            'eval(''kazan(''''AssociateFigure'''',hand, gcf)'');']);
    case {  handles.mToolsShowinFigure,...
            handles.mToolsShowinwindow2,handles.mToolsShowinwindow3,...
            handles.mToolsShowinwindow4,handles.mToolsShowinwindow5, ...
            },

        % Determine the type of plotting
        plot_n   = floor(eventdata / 100 + 0.5);
        vert_dim = floor((eventdata-plot_n*100)/10 + 0.5);
        horz_dim = eventdata-plot_n*100-vert_dim*10;

        % Prepare a figure/axis for the plotting
        if isempty(handles.AssociatedFigure) || h == handles.mToolsShowinFigure
            a = figure;
            set(a,'UserData', struct('MainFigure',handles.MainFigure));
            pMenu = uimenu('Parent', a, 'Label', 'KV_Tools');
            uimenu('Parent', pMenu, 'Label', 'Export Image', 'Callback', 'kv_exportimage(gcf, ''gui'')');  % alsi 12.01.05
            uimenu('Parent', pMenu, 'Label', 'Profile Manager', 'Callback', 'kv_printprofile');  % alsi 12.01.05
            uimenu('Parent', pMenu, 'Label', 'Associate',...
                'Callback', ['ud = get(gcf ,''UserData''); hand = getfield(ud,''MainFigure''); gdata = guidata(hand);',...
                'eval(''kazan(''''AssociateFigure'''',hand, gcf)'');']);

        else
            a = handles.AssociatedFigure;
            set(0, 'CurrentFigure', a);
        end

        switch handles.DataPart
            case {1, 3}
                if vert_dim*horz_dim == 1
                    cla;
                else
                    subplot(vert_dim, horz_dim, plot_n)
                    cla;
                end
                plot(1);
                imh = gca;
            case 2
                clf;
                subplot(2,1,2)
                plot(1);
                imh = gca;
                subplot(2,1,1)
                plot(1);
        end

        %   Plotting
        PlotData(h, gca, imh, a);

        %   workaround for fix scale problem
        axis tight
        axis(imh,'tight')

        % Adding labels
        label = get(handles.tUnits, 'String');
        title(handles.fname, 'Interpreter', 'none')
        ylabel('real');
        xlabel(label);

        if handles.DataPart==2
            subplot(2,1,2);
            ylabel('imag');
            xlabel(label);
        end
        % setting for use a new feature of the ML 7xxx - DataTip
         tver = version;
         MatlabVer = str2num(tver(1:3));
        if MatlabVer >= 7,
             fhandle = @DataTip_Update;
             uimenu('Parent', pMenu, 'Label', 'DataTip X-only',...
                 'Callback', ['dha = datacursormode(gcf); set(dha,''UpdateFcn'',@kv_DataTip_Update);']);
             uimenu('Parent', pMenu, 'Label', 'DataTip Y-only',...
                 'Callback', ['dha = datacursormode(gcf); set(dha,''UpdateFcn'',@kv_DataTip_Update1);']);
             uimenu('Parent', pMenu, 'Label', 'DataTip normal',...
                 'Callback', ['dha = datacursormode(gcf); set(dha,''UpdateFcn'',[]);']);
%              dha = datacursormode(a);
%             set(dha,'UpdateFcn',@kv_DataTip_Update,'SnapToDataVertex','on');
%              datacursormode on
            % mouse click on plot
        end

        %         PaperSize = get(a,'PaperSize');
        %         brd = 0.6;
        %     set(a,'PaperPosition', [brd, brd, PaperSize(1)-2*brd, PaperSize(2)-2*brd]);

    case handles.mToolsSkyline
        h = figure;
        xlim = get(handles.aReal, 'XLim');
        ylim = get(handles.aReal, 'YLim');
        handles.src.ax.Lim = [xlim, ylim] ;
        switch handles.DataPart
            case 1,
                kvskyline(h, handles.src.ax, handles.src.y, 'plot', 'image');
            case 3,
                kvskyline(h, handles.src.ax, abs(handles.src.y), 'plot', 'image');
        end
        uimenu('Parent', h, 'Label', 'Export Image (KV)', 'Callback', 'kv_exportimage(gcf, ''gui'')');  % alsi 12.01.05
    case handles.mToolsSkyline1
        h = figure;
        xlim = get(handles.aReal, 'XLim');
        ylim = get(handles.aReal, 'YLim');
        for ii = GetPlotNumber(handles)
            if ii~=0
                % output
                handles.out{ii}.ax.Lim = [xlim, ylim] ;
                newax = getplotaxis(handles.xrepresent, handles.yrepresent, handles.out{ii}.ax);
                switch handles.DataPart
                    case {1,2}
                        kvskyline(h, newax, handles.out{ii}.y, 'plot', 'contour', 'hold', 'on');
                    case 3,
                        kvskyline(h, newax, abs(handles.out{ii}.y), 'plot', 'contour', 'hold', 'on');
                end
            else
                % source
                handles.src.ax.Lim = [xlim, ylim] ;
                newax = getplotaxis(handles.xrepresent, handles.yrepresent, handles.src.ax);
                switch handles.DataPart
                    case {1,2}
                        kvskyline(h, newax, handles.src.y, 'plot', 'contour');
                    case 3,
                        kvskyline(h, newax, abs(handles.src.y), 'plot', 'contour');
                end
            end
        end
        uimenu('Parent', h, 'Label', 'Export Image (KV)', 'Callback', 'kv_exportimage(gcf, ''gui'')');  % alsi 12.01.05
end
return;


% --------------------------------------------------------------------
function handles = mFileHistory(h, eventdata, handles, varargin)
menus = get(handles.mFile, 'Children');
firstpos = get(handles.mFileNextHistory,'Position')+1;
lastpos = size(menus, 1);
selected = 1;

switch(h)
    case handles.mFileAddHistory
        thisfile = findobj(menus, 'Label', handles.fname);
        if ~isempty(handles.fname) && isempty(thisfile)
            handles.fhistory{end+1}.fname = handles.fname;
            handles.fhistory{end}.ax = handles.src.ax;
        end
        selected = size(handles.fhistory, 2);
    case handles.mFileRemoveHistory
        selmenus = findobj(menus, 'Checked', 'on');
        if isempty(selmenus), return; end;
        delfile = get(selmenus, 'Label');
        for k = 1:size(handles.fhistory, 2)
            if handles.fhistory{k}==delfile
                handles.fhistory = {handles.fhistory{[1:(k-1),(k+1):end]}};
                break;
            end
        end
        selected = 0;
    case handles.mFileRemoveAllHistory
        handles.fhistory = {};
    case handles.mFileNextHistory
        if size(menus, 1) < firstpos, return; end;
        selmenus = findobj(menus, 'Checked', 'on');
        if isempty(selmenus)
            pos = firstpos;
        else
            pos = get(selmenus, 'Position');
            if pos == lastpos, pos = firstpos; else pos = pos+1; end
        end
        nextmenu = findobj(menus, 'Position', pos);
        FileHistory(nextmenu, eventdata, handles, varargin);
        return
    case handles.mFilePreviousHistory
        if size(menus, 1) < firstpos, return;  end;
        selmenus = findobj(menus, 'Checked', 'on');
        if isempty(selmenus)
            pos = firstpos;
        else
            pos = get(selmenus, 'Position');
            if pos == firstpos, pos = lastpos; else pos = pos-1; end
        end
        nextmenu = findobj(menus, 'Position', pos);
        FileHistory(nextmenu, eventdata, handles, varargin);
        return;
end

% menu recreating
delete(findobj(menus, 'UserData', 'history'));

for k=1:size(handles.fhistory, 2)
    if k == 1, sep='on'; else sep = 'off'; end;
    if k == selected, sel='on'; else sel = 'off'; end;
    uimenu(handles.mFile, 'Label',handles.fhistory{k}.fname, 'UserData', 'history',...
        'Tag', 'mPluginsReposition', 'Separator', sep, 'Checked', sel,...
        'Callback', 'kazan(''FileHistory'',gcbo,[],guidata(gcbo))');
end
guidata(h, handles)

% --------------------------------------------------------------------
function FileHistory(h, eventdata, handles, varargin)
menus = get(handles.mFile, 'Children');
set(menus, 'Checked', 'off');
[fpath,name,ext] = fileparts(get(h, 'Label'));
handles.dir_path = fpath;
firstpos = get(handles.mFileNextHistory,'Position')+1;
pos = get(h,'Position') - firstpos + 1;

% set old extension
str = get(handles.pmExt, 'String');
whichext = strcmpi(str, handles.fhistory{pos}.ax.ext);
set(handles.pmExt, 'Value', find(whichext));
handles = pmExt_Callback(gcbo, [], handles, varargin);

str = get(handles.pmDatatype, 'String');
if ~strcmpi(str{get(handles.pmDatatype, 'Value')}, handles.fhistory{pos}.ax.Datatype)
    % set old datatype
    whictype = strcmpi(str, handles.fhistory{pos}.ax.Datatype);
    set(handles.pmDatatype, 'Value', find(whictype));
    handles = pmDatatype_Callback(gcbo, [], handles, varargin);
end

% set old filename
str = get(handles.lDirlist, 'String');
whichfile = strcmp(str, [name, ext]);
set(h, 'Checked', 'on');
set(handles.lDirlist, 'Value', find(whichfile));
lDirlist_Callback(handles.lDirlist, eventdata, handles, varargin)

% --------------------------------------------------------------------
function handles = cmContexMenu1(h, eventdata, handles, varargin)
opt = [handles.cmNormal,handles.cmDiff0, handles.cmDiff1,handles.cmIntegrate];
switch h
    case handles.cmDiff0
        handles.process = 'diff0';
    case handles.cmDiff1
        handles.process = 'diff';
    case handles.cmIntegrate
        handles.process = 'integr';
    otherwise
        handles.process = '';
end
set(opt, 'Checked', 'off');
set(h, 'Checked', 'on');
handles.src.ax.diff = handles.process;
guidata(h,handles);
PlotData(handles.MainFigure);
%lDirlist_Callback(handles.lDirlist, eventdata, handles, varargin);

% --------------------------------------------------------------------
function handles = cmContexMenu2(h, eventdata, handles, varargin)
xopt = [handles.cmxLinear, handles.cmxPoints, handles.cmReverseX, handles.cmSemilogX, handles.cmGspace, ...
    handles.cmxPPM, handles.cmxPPMmin1];
yopt = [handles.cmyLinear, handles.cmyPoints, handles.cmReverseY, handles.cmSemilogY ...
    handles.cmyPPM, handles.cmyPPMmin1];
if h==handles.cmxy
    set(xopt, 'Checked', 'off');
    set(yopt, 'Checked', 'off');
    handles.xrepresent = 'plot';
    handles.yrepresent = 'plot';
    set(handles.cmxLinear, 'Checked', 'on');
    set(handles.cmyLinear, 'Checked', 'on');
elseif sum(h==xopt),
    set(xopt, 'Checked', 'off');
    switch h
        case handles.cmReverseX
            handles.xrepresent = 'reverse';
        case handles.cmSemilogX
            handles.xrepresent = 'semilog';
        case handles.cmGspace
            handles.xrepresent = 'gspace';
        case handles.cmxPPM
            handles.xrepresent = 'ppm';
        case handles.cmxPPMmin1
            handles.xrepresent = 'ppm-1';
        case handles.cmxPoints
            handles.xrepresent = 'pnt';
        otherwise
            handles.xrepresent = 'plot';
    end
    set(h, 'Checked', 'on');
elseif sum(h==yopt),
    set(yopt, 'Checked', 'off');
    switch h
        case handles.cmReverseY
            handles.yrepresent = 'reverse';
        case handles.cmSemilogY
            handles.yrepresent = 'semilog';
        case handles.cmyPPM
            handles.yrepresent = 'ppm';
        case handles.cmyPPMmin1
            handles.yrepresent = 'ppm-1';
        case handles.cmyPoints
            handles.yrepresent = 'pnt';
        otherwise
            handles.represent = 'plot';
    end
    set(h, 'Checked', 'on');
end
guidata(h,handles);
PlotData(handles.MainFigure);

% --------------------------------------------------------------------
function varargout = Filter_CallBack(h, eventdata, handles, varargin)
opt  = [handles.cmNofilter, handles.cmRCfilter, handles.cmMvgaver, handles.cmSavGol];
opt1 = [handles.cmFiltAuto, handles.cmFilt1D, handles.cmFilt2D, handles.cmFilt1D2D];
if sum(h==opt),
    set(opt, 'Checked', 'off');
    switch h
        case handles.cmRCfilter
            handles.filter = 'rc';
        case handles.cmMvgaver
            handles.filter = 'mvgaver';
        case handles.cmSavGol
            handles.filter = 'savgol';
        otherwise
            handles.filter = '';
    end
    set(h, 'Checked', 'on');
    handles.src.ax.filt = handles.filter;
elseif sum(h==opt1),
    set(opt1, 'Checked', 'off');
    switch h
        case handles.cmFilt1D
            handles.filterpar1 = '1D';
        case handles.cmFilt2D
            handles.filterpar1 = '2D';
        case handles.cmFilt1D2D
            handles.filterpar1 = '1D2D';
        otherwise
            handles.filterpar1 = 'auto';
    end
    set(h, 'Checked', 'on');
    handles.src.ax.filtpar1 = handles.filterpar1;
end
guidata(h,handles);
PlotData(handles.MainFigure);
if nargout > 0, varargout{1}=handles; end

% --------------------------------------------------------------------
function mHelp(h, eventdata, handles)

switch h
    case handles.mHelp1,
        msgbox({'Parameters(fields) that affect visualization:'; ...
            'title     - data title'; ...
            'x/y       - 1D axis array (column!)   REQUIRED'; ...
            'x/ylabel  - axis labels'; ...
            'Color     - color of the trace b/g/r/c/m/y/k or [r g b], rgb=0..1';...
            'LineStyle - style of the line -/:/-./--'; ...
            'LineWidth - width of the line 0..inf'; ...
            'Marker    - style of dots none/./o/x/+/*/s/d/v/^/</>/p/h'; ...
            'contour   - for 2D plots - contour level';''; ...
            'Plot type related:'; ...
            'reference - for PPM plot'; ...
            'freq1     - for g-plot';''; ...
            'Data processing:'; ...
            'dx, s          - shift of abscissa/ordinate(2D)'; ...
            'diff, ps_mod    - pseudo modulation parameters diff/integr'; ...
            'filt, tc, pol   - filtering parameters rc/mvgaver/savgol, 1..inf, 1..inf';''; ...
            'Fields of ax structure are accessible throught the';...
            'combobox and edit controls at the bottom of the main window'},...
            'Axis Data Structure (ax)', 'none')
    case handles.mHelp2,
        msgbox({'Press a) middle button or b) left and right'; ...
            'mouse buttons at the same time or ';...
            'c) SHIFT and left mouse button';...
            'to trace the mouse cursor position';'';...
            'In measurements mode:';...
            '2 clicks followed by ENTER - distance';...
            '4 clicks followed by ENTER - ratio'}, 'Mouse Cursor', 'none')
    case handles.mHelp3,
        msgbox({'1: Load data';...
            '    Select file extension';'    Select file type';...
            '    Change loading parameters if needed';...
            '    Find directory';'    Select file';'';...
            '2: Adjust the view';...
            '    Right click on the plot, select options';'';...
            '3: Process the data';...
            '    Run plugin from Plugins menu';...
            '    Do required procesing';'';...
            '4: Copy processed data from Output to Source (<<)';...
            '   Repeat processing if needed';'';...
            '5: Print plot';...
            '   select Tools/Show in Window'},...
            'How To Start', 'none')
    case handles.mHelp4,
        msgbox({'All data-loading routines have similar structure:'; ...
            '[ax,y,dsc] = routine(filename,...),';...
            'where ''ax'' is the axis data structure (see corresponding help),'; ...
            '''y'' is a 2D data array (first dimension corresponds to ''ax.x''),';...
            '''dsc'' is the structure with arbitrary field names (information)'},...
            'Data Load Routines', 'none')
    case handles.mHelp5,
        msgbox({['Plugins are unified data processing routines, which take data ',...
            'from Kazan Viewer and return result to Viewer. The origin of data ',...
            'is called ''Source'' the result is sent to ''Output''. Plugins can be ',...
            'applied in a serial way by transfering of the output of precursed plugin ',...
            'to the Source (button ''<<'' of the Viewer.)'];'';...
            ['Plugins of Kazan Viewer are the "m" modules located in the ',...
            'directory "plugins" of the viewer package. Plugins are loaded to the Viewer',...
            'by a routine "kazan_plugins", to add or remove plugins from the "plugins" menu',...
            'just edit it following the style as such:'];...
            ['  out{end+1}.name = ''Name to appear in the menu''; '];...
            ['  out{end}.file = ''filename_of_the_plugin''; ']; 'WARNING: it is case sensetive';'';...
            ['The Viewer set the ''handles.hh'' field of plugin to the ',...
            'Viewer''s handle . The handle structure of Viewer ',...
            'can then be obtained by ''ha = guidata(handles.hh)'' call.'];'';...
            'The Source data located at ''ha.src''';...
            'The Output data can be stored to ''ha.out{}''';...
            'The ''src'' and ''out'' structures have ''y'', ''ax'' and ''dsc'' fields.';...
            'See Data Load Routines for description.';'';...
            'Function calls:';...
            'The plotting function of plugin have to be finished by:';...
            '    guidata(handles.hh, ha);';...
            '    [path,name,ext,ver] = fileparts(get(h.hh, ''FileName''));';...
            '    eval([name ''(''''SetDataSource'''', 2, hand)'']);';...
            'When data are loaded or plugin is opened';...
            'the ''FileUpdate(handles)'' routine of plugin is called.';'';...
            'Ini file:';...
            'Ini file (kazan.ini) can be accessed using ''inimanaging'' routine.'},...
            'Plugins Programming', 'none')
    case handles.mHelp6,
        msgbox({['Select ''Tools/Show in Separate Figure'' menu.',...
            'Then the standard Matlab figure will appear.'];...
            'Optionally the plot can be conditioned using Export Image utility';...
            'or Profile Manager (see corresponding menus)';'';...
            'To create multiple axis on the figure use Show in Axis menu.';...
            ''},...
            'How to print', 'none')
    otherwise, msgbox('No help found');
end
% --------------------------------------------------------------------
function res = superplot(xrepresent, yrepresent, ax, y, p)
c = safeget(ax, 'Color', 'b');
if ischar(c)
    if strcmp(c, 'jet')||strcmp(c, 'lines')||strcmp(c, 'hsv')
%         ss = {'r', 'g', 'b', 'k', 'm', 'c'};
        cmcc= colormap(c);
        if size(y, 2)<64
            idx = floor(linspace(1, 64, size(y, 2)));
            cmcc = cmcc(idx, :);
        end
        for ii = 1:size(y, 2)
            ncc = floor( 64*(ii/64 - floor(ii/64)) );
            tc{ii} = cmcc(ncc, :);
        end
        c = tc;
    end
end
l = safeget(ax, 'LineStyle', '-');
lw = safeget(ax, 'LineWidth', 1);
m = safeget(ax, 'Marker', 'none');

if iscell(c)
    if length(c)==size(y, 2)
        for ii = 1:length(c)
            res(ii) = plot(ax.x, y(:, ii), 'Color', c{ii}, 'LineStyle', l, 'LineWidth', lw, 'Marker', m, 'Parent', p); %hold on;
        end
    else        
        res = plot(ax.x, y, 'Color', c{1}, 'LineStyle', l, 'LineWidth', lw, 'Marker', m, 'Parent', p);
    end
else
    res = plot(ax.x, y, 'Color', c, 'LineStyle', l, 'LineWidth', lw, 'Marker', m, 'Parent', p);
end

if strcmp(xrepresent, 'semilog'), set(p, 'XScale','log', 'XDir', 'normal');
else
    set(p, 'XScale','linear');
    if strcmp(xrepresent, 'reverse') || strcmp(xrepresent, 'gspace'),
        set(p, 'XDir','reverse');
    else set(p, 'XDir','normal');
    end
end
if strcmp(yrepresent, 'semilog'), set(p, 'YScale','log', 'YDir', 'normal');
else
    set(p, 'YScale','linear');
    if strcmp(yrepresent, 'reverse'), set(p, 'YDir','reverse');
    else set(p, 'YDir','normal');
    end
end
% --------------------------------------------------------------------
function ePar_Callback(h, eventdata, handles, varargin)
str = get(handles.pmSelPar, 'String');
val = get(handles.pmSelPar, 'Value');
fname = str{val}; newstr = get(handles.ePar, 'String');
handles.prop.values{val} = newstr;
if handles.prop.types{val}=='s'
    eval(['sstemp = ''', get(handles.ePar, 'String'), ''';']);
else
    eval(['sstemp = ', get(handles.ePar, 'String'), ';']);
end
handles.src.ax = setfield(handles.src.ax, fname, sstemp);
guidata(handles.MainFigure, handles);
PlotData(handles.MainFigure);

% --------------------------------------------------------------------
function smSelPar_Callback(h, eventdata, handles, varargin)
val = get(handles.pmSelPar, 'Value');
set(handles.ePar, 'String', handles.prop.values{val});

% --------------------------------------------------------------------
function mFile2Source(h, eventdata, handles, varargin)
switch handles.DataSource
    case 1
    case {2, 3}
        switch length(handles.out)
            case 0, n=0;
            case 1, n=1;
            otherwise
                cols = {}; titles = {};
                for k = 1:size(handles.out, 2)
                    cols{end+1, 1} =  safeget(handles.out{k}.ax, 'Color', [0;0;1]);
                    if ischar(cols{end, 1})
                        if cols{end, 1}=='jet'
                            cols{end, 1} = [0, 0, 1];
                        end
                    end
                    titles{end+1, 1} = safeget(handles.out{k}, 'title', '');
                end
                n = kv_seltrace(cols, titles);
        end
        if n
            handles.src = handles.out{n};
            SetDataSource(1, handles);
        end
end

% --------------------------------------------------------------------
function handles = sl2Dplot_Callback(h, eventdata, handles, varargin)
handles.Selection = floor(get(handles.sl2Dplot, 'Value')+0.5);
guidata(h, handles);
str = {'X', 'Y'};
data = safeget(handles.src.ax, lower(str{handles.AntiProjection}), []);
if isempty(data), data = 0; else data = data(min(handles.Selection,length(data))); end
label = safeget(handles.src.ax, [lower(str{handles.AntiProjection}), 'label'], '?');
str1 = [str{handles.AntiProjection}, ' (', label,'): ' num2str(data), ' [', num2str(handles.Selection),']'];
set(handles.sl2Dplot, 'Selected', 'on');
set(handles.stCoordinates, 'String', str1);
PlotData(handles.MainFigure);

% --------------------------------------------------------------------
function pb2DSelect_Callback(h, eventdata, handles, varargin)
str = {'X', 'Y'};
handles.Projection = safeget(handles,'Projection',1);
if handles.Projection == 1, handles.Projection = 2; handles.AntiProjection = 1;
else handles.Projection = 1; handles.AntiProjection = 2;
end
set(h, 'String', str{handles.Projection});
guidata(h, handles);
UpdateSlider(handles);
% --------------------------------------------------------------------
function mFileWorkspace(h, eventdata, handles, varargin)
if h ~= handles.mFileSaveWorksp
    tt = [];
    tt = evalin('base', 'tmp'); %'assignin(''caller'', ''tt'', tmp)');
    if isfield(tt, 'x')
        tt.ax.x = tt.x;
    end
    if ~isfield(tt, 'ax')
        error('tmp structure does not contain "ax" substructure');
    end
    if ~isfield(tt.ax, 'y')
        tt.ax.y = zeros(size(tt.y, 2), 1);
        disp('tmp.ax.y(:) = 0 is added');
    end
    if ~isfield(tt.ax, 'xlabel')
        tt.ax.xlabel = '?, ?';
        disp('tmp.ax.xlabel = ''?, ?'' is added');
    end    
    
    disp('variable ''tmp'' loaded from workspace');
    handles.src = tt;
    guidata(h, handles);
    PlotData(handles.MainFigure);
    if isfield(handles.src.ax, 'title')
        set(handles.MainFigure, 'Name', ['KAZAN - workspace [',handles.src.ax.title,']'])
    else
        set(handles.MainFigure, 'Name', ['KAZAN - workspace [',handles.process,']'])
    end
else
    assignin('base', 'tmp', handles.src);
    disp('variable ''tmp'' created');
    disp(handles.src);
end
% --------------------------------------------------------------------
function Zoom_Callback(h, ev, handles, varargin)
    
lim = get(handles.aReal, 'XLim');
dx  = lim(1)-lim(2);
limy = get(handles.aReal, 'YLim');
dy  = limy(1)-limy(2);

ilim = get(handles.aImag, 'XLim');
idx  = ilim(1)-ilim(2);
ilimy = get(handles.aImag, 'YLim');
idy  = ilimy(1)-ilimy(2);

switch h
    case handles.tbX_half, 
        lim = lim + dx * [1, -1] / 2;
        ilim = ilim + idx * [1, -1] / 2;
    case handles.tbX_twice, 
        lim = lim - dx * [1, -1] / 4;
        ilim = ilim - idx * [1, -1] / 4;
    case handles.tbY_half, 
        limy = limy + dy * [1, -1] / 2;
        ilimy = ilimy + idy * [1, -1] / 2;
    case handles.tbY_twice, 
        limy = limy - dy * [1, -1] / 4;
        ilimy = ilimy - idy * [1, -1] / 4;        
    case handles.tbPan, 
        if get(handles.tbPan, 'Value')
            pan(handles.MainFigure, 'on');
        else
            pan(handles.MainFigure, 'off');
        end
    case handles.tbImReX
        ilim = lim;
end
set(handles.aReal, 'XLim', lim, 'YLim', limy);
set(handles.aImag, 'XLim', ilim, 'YLim', ilimy );

% --------------------------------------------------------------------
function mViewRange(h, eventdata, handles, varargin)
viewarea_func(handles)

% --------------------------------------------------------------------
function newax = getplotaxis(xfname, yfname, ax)
newax = ax;

dx = safeget(ax, 'dx', 0); newax.dx=0;
dy = safeget(ax, 'dy', 0); newax.dy=0;
gref = safeget(ax, 'freq1', 34E9)*planck/bmagn*1e4;
if gref == 0.0, gref = 1.; end
ref = safeget(ax, 'reference', 34E9);
if ref == 0.0, ref = 1.; end

xlabel = safeget(ax, 'xlabel', '');
pos = findstr(xlabel, ',');
if isempty(pos), xvar='?,'; else xvar = xlabel(1:min(pos)); end
ylabel = safeget(ax, 'ylabel', '');
pos = findstr(ylabel, ',');
if isempty(pos), yvar='?,'; else yvar = ylabel(1:min(pos)); end

if strcmp(xfname, 'gspace')
    newax.x = gref./(ax.x+dx);
    newax.xlabel = 'g-factor,';
elseif strcmp(xfname, 'ppm-1')
    newax.x = ((ax.x+dx)./ref - 1)*1E6;
    newax.xlabel = [xlabel,'ppm'];
elseif strcmp(xfname, 'ppm')
    newax.x = (ax.x+dx)./ref*1E6;
    newax.xlabel = [xvar,'ppm'];
elseif strcmp(xfname, 'pnt')
    newax.x = (1:size(ax.x,1)).';
    newax.xlabel = [xvar,'pts'];
else
    newax.x = ax.x+dx;
    newax.xlabel = xlabel;
end

if strcmp(yfname, 'gspace')
    newax.y = gref./(ax.y+dy);
    newax.ylabel = 'g-factor,';
elseif strcmp(yfname, 'ppm-1')
    newax.y = ((ax.y+dy)./ref - 1)*1E6;
    newax.ylabel = [yvar,'ppm'];
elseif strcmp(yfname, 'ppm')
    newax.y = (ax.y+dy)./ref*1E6;
    newax.ylabel = [yvar,'ppm'];
elseif strcmp(yfname, 'pnt')
    newax.y = (1:size(ax.y,1))';
    newax.ylabel = [yvar,'pts'];
else
    newax.y = ax.y+dy;
    newax.ylabel = ylabel;
end

% ---------------------------------------------------
function varargout = getcoord(varargin)
% Calculates cursor position in axes units using 'CurrentPoint' of
% CurrentAxes. if nargout=3, return handle of this axes
handles = varargin{1};
hand = [];
if nargin == 2, hand = varargin{2}; end
% Determination of positon of the pointer
CoordXY = get(handles.MainFigure, 'CurrentPoint');
RePos = get(handles.aReal, 'Position');
ImPos = get(handles.aImag, 'Position');
ReX = 0;
ReY = 0;
ImX = 0;
ImY = 0;
if isempty(hand)
    switch handles.DataPart
        case {2}
            if CoordXY(1)>RePos(1)&& CoordXY(1)< sum(RePos([1, 3])), ReX = 1; end
            if CoordXY(2)>RePos(2)&& CoordXY(2)< sum(RePos([2, 4])), ReY = 1; end
            if CoordXY(1)>ImPos(1)&& CoordXY(1)< sum(ImPos([1, 3])), ImX = 1; end
            if CoordXY(2)>ImPos(2)&& CoordXY(2)< sum(ImPos([2, 4])), ImY = 1; end
        case {1,3}
            if CoordXY(1)>RePos(1)&& CoordXY(1)< sum(RePos([1, 3])), ReX = 1; end
            if CoordXY(2)>RePos(2)&& CoordXY(2)< sum(RePos([2, 4])), ReY = 1; end
    end
    if isempty(hand),
        if ReX&&ReY, hand = handles.aReal; end
        if ImX&&ImY, hand = handles.aImag; end
    end
end
x=[]; y=[];
if ~isempty(hand),
    CXY = get(hand, 'CurrentPoint');
    x = CXY(1, 1);
    y = CXY(1, 2);
end
varargout{1} = x;
varargout{2} = y;
if nargout == 3,
    varargout{3} = hand;
end

% ---------------------------------------------------
function rs = GetSlice(src, dim, sel)
src.ax.y = safeget(src.ax,'y',0);
rs = src;
switch dim
    case 1 % X
        rs.y    = src.y(:,min(sel,end));
        rs.ax.y = src.ax.y(min(sel,end));
    case 2 % Y
        rs.y    = src.y(sel,:)';
        rs.ax.x = src.ax.y;
        rs.ax.y = src.ax.x(min(sel,end));
        rs.ax.xlabel = safeget(src.ax,'ylabel', '?,');
        rs.ax.ylabel = safeget(src.ax,'xlabel', '?,');
    otherwise
end

% ---------------------------------------------------
function rs = GetSourceSlice(handles)
rs = GetSlice(handles.src, handles.Projection, handles.Selection);

% ---------------------------------------------------
function UpdateSlider(handles)
n = 1; % number of points
switch handles.DataSource
    case 1
        n = size(handles.src.y, handles.AntiProjection);
    case 2
        for jj=1:length(handles.out)
            n1 = size(handles.out{jj}.y, handles.AntiProjection);
            if n1>n, n=n1; end
        end
    otherwise
        n = size(handles.src.y, handles.AntiProjection);
        for jj=1:length(handles.out)
            n1 = size(handles.out{jj}.y, handles.AntiProjection);
            if n1>n, n=n1; end
        end
end
if n > 1
    nold = get(handles.sl2Dplot, 'Max');
    if nold ~= n
        set(handles.sl2Dplot, 'Min', 1, 'Max', n, 'Value', 1, 'SliderStep', [1,1]/(n-1), 'Enable', 'on');
    else
        set(handles.sl2Dplot, 'Enable', 'on');
    end
else
    set(handles.sl2Dplot, 'Min', 1, 'Max', 2, 'Value', 1, 'SliderStep', [1,1], 'Enable', 'off');
end
sl2Dplot_Callback(handles.sl2Dplot, [], handles);

% ---------------------------------------------------
function [x,y, button]=GetPoint(handles, n)
set(0,'CurrentFigure', handles.MainFigure)
[x,y, button]=ginput(n);

% ---------------------------------------------------
function [x,y]=GetBox(handles)
set(0,'CurrentFigure', handles.MainFigure)
waitforbuttonpress;
[x1, y1, hand] = getcoord(handles);
rbbox;
[x2, y2] = getcoord(handles,hand);
x = [x1; x2];
y = [y1; y2];

% ---------------------------------------------------
function fMoveDataBr_Callback(h, handles)
handles.changePart = 1;
guidata(h, handles);

% ---------------------------------------------------
function AssociateFigure(h, a_fig)
hand = guidata(h);
hand.AssociatedFigure = a_fig;
if isempty(a_fig)
    set(hand.mToolsAssociatedwindow, 'Label', 'No Figure Is Associated');
else
    set(hand.mToolsAssociatedwindow, 'Label', ['Associated Figure -> ',num2str(gcf)]);
end
guidata(h,hand);


