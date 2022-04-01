handles.mLoadOneByOne = uimenu(handles.figure1, 'Tag', 'mLoadOneByOne', 'Label', 'Load Data', 'Checked', 'off', 'Separator', 'off', 'Callback', 'fittingbox('butLoadData_Callback',gcbo,[],guidata(gcbo))' ); 
handles. = uimenu(handles.figure1, 'Tag', '', 'Label', 'Simulation', 'Checked', 'off', 'Separator', 'off', 'Callback', '' ); 
handles.mSimFFunc = uimenu(handles., 'Tag', 'mSimFFunc', 'Label', 'Select Function ...', 'Checked', 'off', 'Separator', 'off', 'Callback', 'fittingbox('mSimFFunc_Callback',gcbo,[],guidata(gcbo))', 'Accelerator', '' ); 
handles.mSimAgr = uimenu(handles., 'Tag', 'mSimAgr', 'Label', 'Source processing', 'Checked', 'on', 'Separator', 'on', 'Callback', 'fittingbox('mSimAgr_Callback',gcbo,[],guidata(gcbo))', 'Accelerator', '' ); 
handles.mSimDif = uimenu(handles., 'Tag', 'mSimDif', 'Label', 'Show Difference', 'Checked', 'on', 'Separator', 'off', 'Callback', 'fittingbox('mSimDif_Callback',gcbo,[],guidata(gcbo))', 'Accelerator', '' ); 
handles.mSimShowBds = uimenu(handles., 'Tag', 'mSimShowBds', 'Label', 'Show Range Bounds', 'Checked', 'off', 'Separator', 'off', 'Callback', 'fittingbox('mSimShowBds_Callback',gcbo,[],guidata(gcbo))', 'Accelerator', '' ); 
handles.mSimAddRng = uimenu(handles., 'Tag', 'mSimAddRng', 'Label', 'Add Simulation Range', 'Checked', 'off', 'Separator', 'on', 'Callback', 'fittingbox('mSimAddRng_Callback',gcbo,[],guidata(gcbo))', 'Accelerator', '' ); 
handles.mSimDelRng = uimenu(handles., 'Tag', 'mSimDelRng', 'Label', 'Delete Simulation Range', 'Checked', 'off', 'Separator', 'off', 'Callback', 'fittingbox('mSimDelRng_Callback',gcbo,[],guidata(gcbo))', 'Accelerator', '' ); 
handles.edNIter = uicontrol(handles.figure1, 'Style', 'edit', 'Tag', 'edNIter', 'String', '100', 'Value', 0, 'Units', 'pixels', 'Callback', 'fittingbox('edNIter_Callback',gcbo,[],guidata(gcbo))'); 
handles.txtNIter = uicontrol(handles.figure1, 'Style', 'text', 'Tag', 'txtNIter', 'String', 'N. iterations', 'Value', 0, 'Units', 'pixels', 'Callback', ''); 
handles.pbShowHide = uicontrol(handles.figure1, 'Style', 'pushbutton', 'Tag', 'pbShowHide', 'String', 'v', 'Value', 0, 'Units', 'pixels', 'Callback', 'fittingbox('pbShoHide_Callback',gcbo,[],guidata(gcbo))'); 
handles.edRangeVal = uicontrol(handles.figure1, 'Style', 'edit', 'Tag', 'edRangeVal', 'String', '0', 'Value', 0, 'Units', 'pixels', 'Callback', 'fittingbox('edRangeVal_Callback',gcbo,[],guidata(gcbo))'); 
handles.slRangeStep = uicontrol(handles.figure1, 'Style', 'slider', 'Tag', 'slRangeStep', 'String', '{, ''}', 'Value', 0, 'Units', 'pixels', 'Callback', 'fittingbox('slRangeStep_Callback',gcbo,[],guidata(gcbo))'); 
handles.txtDatas = uicontrol(handles.figure1, 'Style', 'text', 'Tag', 'txtDatas', 'String', 'Data selector', 'Value', 0, 'Units', 'pixels', 'Callback', ''); 
handles.txtMain = uicontrol(handles.figure1, 'Style', 'text', 'Tag', 'txtMain', 'String', 'Fitting settings', 'Value', 0, 'Units', 'pixels', 'Callback', ''); 
handles.txtRange = uicontrol(handles.figure1, 'Style', 'text', 'Tag', 'txtRange', 'String', 'Fitting range', 'Value', 0, 'Units', 'pixels', 'Callback', ''); 
handles.pbKillRange = uicontrol(handles.figure1, 'Style', 'pushbutton', 'Tag', 'pbKillRange', 'String', '-', 'Value', 0, 'Units', 'pixels', 'Callback', 'fittingbox('mSimDelRng_Callback',gcbo,[],guidata(gcbo))'); 
handles.pbAddRange = uicontrol(handles.figure1, 'Style', 'pushbutton', 'Tag', 'pbAddRange', 'String', '+', 'Value', 0, 'Units', 'pixels', 'Callback', 'fittingbox('mSimAddRng_Callback',gcbo,[],guidata(gcbo))'); 
handles.listbox2 = uicontrol(handles.figure1, 'Style', 'listbox', 'Tag', 'listbox2', 'String', '', 'Value', 1, 'Units', 'pixels', 'Callback', 'fittingbox('listbox2_Callback',gcbo,[],guidata(gcbo))'); 
handles.frRange = uicontrol(handles.figure1, 'Style', 'frame', 'Tag', 'frRange', 'String', '', 'Value', 0, 'Units', 'pixels', 'Callback', '0;'); 
handles.popSel = uicontrol(handles.figure1, 'Style', 'popupmenu', 'Tag', 'popSel', 'String', '1', 'Value', 1, 'Units', 'pixels', 'Callback', 'fittingbox('popSel_Callback',gcbo,[],guidata(gcbo))'); 
handles.butFittAll = uicontrol(handles.figure1, 'Style', 'pushbutton', 'Tag', 'butFittAll', 'String', 'Fitt all', 'Value', 0, 'Units', 'pixels', 'Callback', 'fittingbox('butFittAll_Callback',gcbo,[],guidata(gcbo))'); 
handles.but1Step = uicontrol(handles.figure1, 'Style', 'pushbutton', 'Tag', 'but1Step', 'String', '1 step', 'Value', 0, 'Units', 'pixels', 'Callback', 'fittingbox('but1Step_Callback',gcbo,[],guidata(gcbo))'); 
handles.butFitt = uicontrol(handles.figure1, 'Style', 'pushbutton', 'Tag', 'butFitt', 'String', 'Fitt', 'Value', 0, 'Units', 'pixels', 'Callback', 'fittingbox('butFitt_Callback',gcbo,[],guidata(gcbo))'); 
handles.butCol = uicontrol(handles.figure1, 'Style', 'pushbutton', 'Tag', 'butCol', 'String', '', 'Value', 0, 'Units', 'pixels', 'Callback', 'fittingbox('butCol_Callback',gcbo,[],guidata(gcbo))'); 
handles.popup = uicontrol(handles.figure1, 'Style', 'popupmenu', 'Tag', 'popup', 'String', '{, '1e-6', '1e-5', '1e-4', '1e-3', '1e-2', '0.1', '1', '10', '100', '1e+3', '1e+4', '1e+5', '1e+6'}', 'Value', 6, 'Units', 'pixels', 'Callback', 'fittingbox('popup_Callback',gcbo,[],guidata(gcbo))'); 
handles.slPars = uicontrol(handles.figure1, 'Style', 'slider', 'Tag', 'slPars', 'String', '{, ''}', 'Value', 0, 'Units', 'pixels', 'Callback', 'fittingbox('slPars_Callback',gcbo,[],guidata(gcbo))'); 
handles.ePars = uicontrol(handles.figure1, 'Style', 'edit', 'Tag', 'ePars', 'String', 'polyfit(x, y, n) ', 'Value', 0, 'Units', 'pixels', 'Callback', 'fittingbox('ePars_Callback',gcbo,[],guidata(gcbo))'); 
handles.list = uicontrol(handles.figure1, 'Style', 'listbox', 'Tag', 'list', 'String', '{, 'f = polyfit(x, y, n)', 'n = 2'}', 'Value', 1, 'Units', 'pixels', 'Callback', 'fittingbox('list_Callback',gcbo,[],guidata(gcbo))'); 
handles.frMain = uicontrol(handles.figure1, 'Style', 'frame', 'Tag', 'frMain', 'String', '{, ''}', 'Value', 0, 'Units', 'pixels', 'Callback', ''); 
handles.frChoose = uicontrol(handles.figure1, 'Style', 'frame', 'Tag', 'frChoose', 'String', '{, ''}', 'Value', 0, 'Units', 'pixels', 'Callback', ''); 
