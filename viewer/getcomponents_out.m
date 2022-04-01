handles.Untitled_1 = uimenu(handles.figure1, 'Tag', 'Untitled_1', 'Label', 'About...', 'Checked', 'off', 'Separator', 'off', 'Callback', '' ); 
handles.mAboutHelp = uimenu(handles.Untitled_1, 'Tag', 'mAboutHelp', 'Label', 'Help', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('mAbout_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mAbout = uimenu(handles.Untitled_1, 'Tag', 'mAbout', 'Label', 'About', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('mAbout_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mOptions = uimenu(handles.figure1, 'Tag', 'mOptions', 'Label', 'Options', 'Checked', 'off', 'Separator', 'off', 'Callback', '' ); 
handles. = uimenu(handles.mOptions, 'Tag', '', 'Label', 'Isotopes (EasySpin)', 'Checked', 'off', 'Separator', 'off', 'Callback', 'isotopes' ); 
handles.mSimLoadRange = uimenu(handles.mOptions, 'Tag', 'mSimLoadRange', 'Label', 'Load Range', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('mSimLoadRange_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mSimApplyFreq = uimenu(handles.mOptions, 'Tag', 'mSimApplyFreq', 'Label', 'Auto Load Frequency', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('mSimApplyFreq_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mSimLoadFreq = uimenu(handles.mOptions, 'Tag', 'mSimLoadFreq', 'Label', 'Load Frequency', 'Checked', 'off', 'Separator', 'on', 'Callback', 'simplugin('mSimLoadFreq_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mOptAddScripts = uimenu(handles.mOptions, 'Tag', 'mOptAddScripts', 'Label', 'Add Scripts to DSC', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('mCheckOffOn_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mSimAgr = uimenu(handles.mOptions, 'Tag', 'mSimAgr', 'Label', 'Apply Source Processing', 'Checked', 'on', 'Separator', 'off', 'Callback', 'simplugin('mSimAgr_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mCalc = uimenu(handles.mOptions, 'Tag', 'mCalc', 'Label', 'Calculate', 'Checked', 'off', 'Separator', 'on', 'Callback', 'simplugin('check_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mSimLarm = uimenu(handles.mOptions, 'Tag', 'mSimLarm', 'Label', 'Shift to Larmor Frequency', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('mSimLarm_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mSimChi = uimenu(handles.mOptions, 'Tag', 'mSimChi', 'Label', 'Show Chi^2', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('mSimChi_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mColDiff = uimenu(handles.mOptions, 'Tag', 'mColDiff', 'Label', '"Diff" color', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('mSimCol_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mColSum = uimenu(handles.mOptions, 'Tag', 'mColSum', 'Label', '"Sum" color', 'Checked', 'off', 'Separator', 'on', 'Callback', 'simplugin('mSimCol_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mSimAutoBaseline = uimenu(handles.mOptions, 'Tag', 'mSimAutoBaseline', 'Label', 'Auto baseline', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('ZZ_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mSimSubBl = uimenu(handles.mOptions, 'Tag', 'mSimSubBl', 'Label', 'Subtract baseline', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('mSimShow_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mShdiff = uimenu(handles.mOptions, 'Tag', 'mShdiff', 'Label', 'Show diff', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('mSimShow_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mShsumm = uimenu(handles.mOptions, 'Tag', 'mShsumm', 'Label', 'Show sum', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('mSimShow_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mShsim = uimenu(handles.mOptions, 'Tag', 'mShsim', 'Label', 'Show sim', 'Checked', 'on', 'Separator', 'off', 'Callback', 'simplugin('mShsim_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mContext1 = uicontextmenu(handles.figure1, 'Tag', 'mContext1', 'Callback', 'simplugin('mContext1_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mCont1AddRange = uicontextmenu(handles.mContext1, 'Tag', 'mCont1AddRange', 'Label', 'Add Range', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('mCont1AddRange_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mContext = uicontextmenu(handles.figure1, 'Tag', 'mContext', 'Callback', '' ); 
handles.mMoveDown = uicontextmenu(handles.mContext, 'Tag', 'mMoveDown', 'Label', 'Move Down', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('mFileMove_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mMoveUp = uicontextmenu(handles.mContext, 'Tag', 'mMoveUp', 'Label', 'MoveUp', 'Checked', 'off', 'Separator', 'on', 'Callback', 'simplugin('mFileMove_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mConMakeLoc = uicontextmenu(handles.mContext, 'Tag', 'mConMakeLoc', 'Label', 'Make it local', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('mConComm_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mConMakeCom = uicontextmenu(handles.mContext, 'Tag', 'mConMakeCom', 'Label', 'Make it common', 'Checked', 'off', 'Separator', 'on', 'Callback', 'simplugin('mConComm_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mContDelPar = uicontextmenu(handles.mContext, 'Tag', 'mContDelPar', 'Label', 'Change Parameter', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('mFileChPar_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mContChPar = uicontextmenu(handles.mContext, 'Tag', 'mContChPar', 'Label', 'Delete Parameter', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('mFileDelPar_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mContAddPar = uicontextmenu(handles.mContext, 'Tag', 'mContAddPar', 'Label', 'Add Parameters', 'Checked', 'off', 'Separator', 'on', 'Callback', 'simplugin('mFileAddPar_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mContNotOpt = uicontextmenu(handles.mContext, 'Tag', 'mContNotOpt', 'Label', 'Do not optimise', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('mContOptSet_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mContOpt = uicontextmenu(handles.mContext, 'Tag', 'mContOpt', 'Label', 'Optimise', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('mContOptSet_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mEdit = uimenu(handles.figure1, 'Tag', 'mEdit', 'Label', 'Edit', 'Checked', 'off', 'Separator', 'off', 'Callback', '' ); 
handles.mRefresh = uimenu(handles.mEdit, 'Tag', 'mRefresh', 'Label', 'Refresh', 'Checked', 'off', 'Separator', 'on', 'Callback', 'simplugin('mRefresh_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mFileChName = uimenu(handles.mEdit, 'Tag', 'mFileChName', 'Label', 'Simulation Name ...', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('mFileChName_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mSortPars = uimenu(handles.mEdit, 'Tag', 'mSortPars', 'Label', 'Sort Parameters', 'Checked', 'off', 'Separator', 'on', 'Callback', 'simplugin('mSortPars_Callback',gcbo,[],guidata(gcbo))' ); 
handles. = uimenu(handles.mEdit, 'Tag', '', 'Label', 'Script Modification', 'Checked', 'off', 'Separator', 'on', 'Callback', '' ); 
handles.mFileOpenEditor = uimenu(handles.mEdit, 'Tag', 'mFileOpenEditor', 'Label', 'Open in Editor', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('mFileOpenEditor_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mFileCopyClipboard = uimenu(handles.mEdit, 'Tag', 'mFileCopyClipboard', 'Label', 'Copy to Clipboard', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('mSave_Callback',gcbo,[],guidata(gcbo))' ); 
handles. = uimenu(handles.figure1, 'Tag', '', 'Label', 'File', 'Checked', 'off', 'Separator', 'off', 'Callback', '' ); 
handles.mDeleteAll = uimenu(handles., 'Tag', 'mDeleteAll', 'Label', 'Delete All', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('butDel_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mFileDelete = uimenu(handles., 'Tag', 'mFileDelete', 'Label', 'Delete', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('butDel_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mSaveAll = uimenu(handles., 'Tag', 'mSaveAll', 'Label', 'Save All ...', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('mSave_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mSave = uimenu(handles., 'Tag', 'mSave', 'Label', 'Save ...', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('mSave_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mFileReload = uimenu(handles., 'Tag', 'mFileReload', 'Label', 'Reload', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('mFileReload_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mFileCopy = uimenu(handles., 'Tag', 'mFileCopy', 'Label', 'Copy', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('mFileCopy_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mLoadAll = uimenu(handles., 'Tag', 'mLoadAll', 'Label', 'Load All...', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('mLoad_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mLoad = uimenu(handles., 'Tag', 'mLoad', 'Label', 'Load ...', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('mLoad_Callback',gcbo,[],guidata(gcbo))' ); 
handles.mFileBaseline = uimenu(handles., 'Tag', 'mFileBaseline', 'Label', 'Add baseline ...', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('AddSim_Callback',gcbo,[],guidata(gcbo))' ); 
handles.AddSim = uimenu(handles., 'Tag', 'AddSim', 'Label', 'Add sim ...', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin('AddSim_Callback',gcbo,[],guidata(gcbo))' ); 
handles.pbReCalc = uicontrol(handles.figure1, 'Style', 'pushbutton', 'Tag', 'pbReCalc', 'String', '', 'Value', 0, 'Units', 'pixels', 'Callback', 'simplugin('pbReCalc_Callback',gcbo,[],guidata(gcbo))'); 
handles.pbCalcAll = uicontrol(handles.figure1, 'Style', 'pushbutton', 'Tag', 'pbCalcAll', 'String', '', 'Value', 0, 'Units', 'pixels', 'Callback', 'simplugin('calc_all', guidata(gcbo))'); 
handles.list = uicontrol(handles.figure1, 'Style', 'listbox', 'Tag', 'list', 'String', '{, 'l-Sys.S = 0.5', 'l-Sys.g(1) = 2', 'l-Sys.g(2) = 2', 'l-Sys.g(3) = 2', 'l-Sys.lw = 5', 'l-Exp.Range = [330,370]', 'l-Sys.gStrain(1) = 0', 'l-Sys.gStrain(2) = 0', 'l-Sys.gStrain(3) = 0', 'l-Exp.mwFreq = 9.77', 'l-Exp.Harmonic = 1', 'l-Shift.x = 0', 'l-Shift.Scale = 1'}', 'Value', 1, 'Units', 'pixels', 'Callback', 'simplugin('list_Callback',gcbo,[],guidata(gcbo))'); 
handles.text4 = uicontrol(handles.figure1, 'Style', 'text', 'Tag', 'text4', 'String', 'Simulation', 'Value', 0, 'Units', 'pixels', 'Callback', ''); 
handles.lbOpt = uicontrol(handles.figure1, 'Style', 'listbox', 'Tag', 'lbOpt', 'String', '{, ' -Opm.Amp = 1', ' -Opm.wf = 0', ' -Opm.dys = 0'}', 'Value', 1, 'Units', 'pixels', 'Callback', 'simplugin('lbOpt_Callback',gcbo,[],guidata(gcbo))'); 
handles.pmOptMethod = uicontrol(handles.figure1, 'Style', 'popupmenu', 'Tag', 'pmOptMethod', 'String', '{, 'None', 'Fit Peak-To-Peak', 'Fit Amplitude', 'Fit Integral', 'Fit Double Integral', 'Fix Double Integral', 'Default optimization', '10 Iterations', '100 Iterations', '1.0E-6 tolerance', '1.0E-8 tolerance'}', 'Value', 2, 'Units', 'pixels', 'Callback', 'simplugin('ZZ_Callback',gcbo,[],guidata(gcbo))'); 
handles.butCalc = uicontrol(handles.figure1, 'Style', 'togglebutton', 'Tag', 'butCalc', 'String', '', 'Value', 0, 'Units', 'pixels', 'Callback', 'simplugin('butCalc_Callback',gcbo,[],guidata(gcbo))'); 
handles.butCol = uicontrol(handles.figure1, 'Style', 'pushbutton', 'Tag', 'butCol', 'String', '', 'Value', 0, 'Units', 'pixels', 'Callback', 'simplugin('mSimCol_Callback',gcbo,[],guidata(gcbo))'); 
handles.popSel = uicontrol(handles.figure1, 'Style', 'popupmenu', 'Tag', 'popSel', 'String', '1', 'Value', 1, 'Units', 'pixels', 'Callback', 'simplugin('popSel_Callback',gcbo,[],guidata(gcbo))'); 
handles.eOpt = uicontrol(handles.figure1, 'Style', 'edit', 'Tag', 'eOpt', 'String', '0', 'Value', 0, 'Units', 'pixels', 'Callback', 'simplugin('eOpt_Callback',gcbo,[],guidata(gcbo))'); 
handles.pmOpt = uicontrol(handles.figure1, 'Style', 'popupmenu', 'Tag', 'pmOpt', 'String', '{, '1e-6', '1e-5', '1e-4', '1e-3', '1e-2', '0.1', '1', '10', '1.0e2', '1.0e3', '1.0e4', '1.0e5', '1.0e6', '1.0e7', '1.0e8'}', 'Value', 6, 'Units', 'pixels', 'Callback', 'simplugin('pmOpt_Callback',gcbo,[],guidata(gcbo))'); 
handles.slOpt = uicontrol(handles.figure1, 'Style', 'slider', 'Tag', 'slOpt', 'String', '{, ''}', 'Value', 0, 'Units', 'pixels', 'Callback', 'simplugin('slOpt_Callback',gcbo,[],guidata(gcbo))'); 
handles.popup = uicontrol(handles.figure1, 'Style', 'popupmenu', 'Tag', 'popup', 'String', '{, '1e-6', '1e-5', '1e-4', '1e-3', '1e-2', '0.1', '1', '10', '100', '1e+3', '1e+4', '1e+5', '1e+6'}', 'Value', 6, 'Units', 'pixels', 'Callback', 'simplugin('popup_Callback',gcbo,[],guidata(gcbo))'); 
handles.slPars = uicontrol(handles.figure1, 'Style', 'slider', 'Tag', 'slPars', 'String', '{, ''}', 'Value', 0, 'Units', 'pixels', 'Callback', 'simplugin('slPars_Callback',gcbo,[],guidata(gcbo))'); 
handles.ePars = uicontrol(handles.figure1, 'Style', 'edit', 'Tag', 'ePars', 'String', '0', 'Value', 0, 'Units', 'pixels', 'Callback', 'simplugin('ePars_Callback',gcbo,[],guidata(gcbo))'); 
handles.frame1 = uicontrol(handles.figure1, 'Style', 'frame', 'Tag', 'frame1', 'String', '{, ''}', 'Value', 0, 'Units', 'pixels', 'Callback', ''); 
handles.text3 = uicontrol(handles.figure1, 'Style', 'text', 'Tag', 'text3', 'String', 'Optimization', 'Value', 0, 'Units', 'pixels', 'Callback', ''); 
handles.frame2 = uicontrol(handles.figure1, 'Style', 'frame', 'Tag', 'frame2', 'String', '{, ''}', 'Value', 0, 'Units', 'pixels', 'Callback', ''); 
