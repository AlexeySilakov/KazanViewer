function varargout = fftfilter(varargin)
% FFTFILTER Application M-file for fftfilter.fig
%    FIG = FFTFILTER launch fftfilter GUI.
%    FFTFILTER('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 22-Jan-2006 12:32:28

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'new');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
    handles.hh = 0;
    handles.Band = [];
    handles = SetMode(handles, 0);
    handles = UpdateList(handles);
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


% --------------------------------------------------------------------
function FileUpdate(varargin);
return;
% -------------------------------------------------------------------
function ret = OnOff(arg)
if arg < 0.5, ret = 'off';
else ret = 'on';
end
   
% -------------------------------------------------------------------
function handles = SetMode(handles, mode)
handles.Mode = mode;

set(handles.rbData,'Enable', OnOff(mode~=0))
set(handles.rbFrequency,'Enable', OnOff(mode~=0))
set(handles.pbAddRange,'Enable', OnOff(mode==2))

set(handles.rbData,'Value', mode==1)
set(handles.rbFrequency,'Value', mode==2)

guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function varargout = lbParsList_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = slPars_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = pmBandSelector_Callback(h, eventdata, handles, varargin)
handles = UpdateList(handles);

% --------------------------------------------------------------------
function varargout = edit1_Callback(h, eventdata, handles, varargin)


% --------------------------------------------------------------------
function varargout = pbApply_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = pbAddRange_Callback(h, eventdata, handles, varargin)
hand = guidata(handles.hh);

axx =hand.out{1}.ax.x;
axy =hand.out{1}.y;

[path,name,ext] = fileparts(get(handles.hh, 'FileName'));
[x,y]=eval([name '(''GetPoint'', hand, 2)']);
if length(x)~=2 return; end
    
[aa,pos1] = min(abs(axx - x(1)));
[aa,pos2] = min(abs(axx - x(2)));

handles.Band(end+1).Pos1 = pos1;
handles.Band(end).Pos2 = pos2;

handles = UpdateList(handles);

Bands   = get(handles.pmBandSelector, 'String');
set(handles.pmBandSelector, 'Value', length(Bands));

guidata(handles.figure1, handles);
Plot(handles);

% --------------------------------------------------------------------
function varargout = pbDeleteRange_Callback(h, eventdata, handles, varargin)
BandIdx = get(handles.pmBandSelector, 'Value');

Band = [];
for ii=1:length(handles.Band)
    if ii~=BandIdx
        if isempty(Band)  
            Band = handles.Band(ii);
        else
            Band(end+1) = handles.Band(ii);
        end
    end
end
handles.Band = Band;

handles = UpdateList(handles);
guidata(handles.figure1, handles);
Plot(handles);


% --------------------------------------------------------------------
function varargout = rbFrequency_Callback(h, eventdata, handles, varargin)

handles = SetMode(handles, 2);
Plot(handles);


% --------------------------------------------------------------------
function varargout = rbData_Callback(h, eventdata, handles, varargin)

handles = SetMode(handles, 1);
Plot(handles);


% --------------------------------------------------------------------
function varargout = pbSetData_Callback(h, eventdata, handles, varargin)
hand = guidata(handles.hh);
handles.Data = hand.src;
handles = SetMode(handles, 1);
guidata(handles.figure1, handles);

Plot(handles);

% --------------------------------------------------------------------
function handles = Plot(handles)
handles = guidata(handles.figure1);
hand = guidata(handles.hh);
hand.out = [];

[sz1, sz2] = size(handles.Data.y);

switch(handles.Mode)
case 0,
case 1,
    % Data
    hand.out{1}.ax = handles.Data.ax;
    hand.out{1}.ax.Color = [0 0 1];
    hand.out{1}.ax.s = 0;
    hand.out{1}.y = handles.Data.y;
    hand.out{1}.title = 'Src';
    
    yfft  = fft(handles.Data.y);
    filf_fft = yfft;
    for ii=1:length(handles.Band)
        filf_fft([handles.Band(ii).Pos1:handles.Band(ii).Pos2],:) = 0;
        filf_fft(sz1-[handles.Band(ii).Pos1:handles.Band(ii).Pos2],:) = 0;
    end
    
    y1 = ifft(filf_fft);

    % Filtered spectrum
    hand.out{2}.ax = handles.Data.ax;
    hand.out{2}.ax.Color = [1 0 0];
    hand.out{2}.ax.s = 0;
    hand.out{2}.y = y1;
    hand.out{2}.title = 'Filtered';
case 2,
    yfft  = fft(handles.Data.y);
    fftx = [1:sz1/2];
    % Data
    hand.out{1}.ax = handles.Data.ax;
    hand.out{1}.ax.Color = [0 0 1];
    hand.out{1}.ax.s = 0;
    hand.out{1}.ax.x = fftx;
    hand.out{1}.y = abs(yfft(1:sz1/2,:));
    hand.out{1}.title = 'Src';
    
    filf_fft = yfft;
    for ii=1:length(handles.Band)
        filf_fft([handles.Band(ii).Pos1:handles.Band(ii).Pos2],:) = 0;
        filf_fft(sz1-[handles.Band(ii).Pos1:handles.Band(ii).Pos2],:) = 0;
    end
    
    % Filtered spectrum
    hand.out{2}.ax = handles.Data.ax;
    hand.out{2}.ax.Color = [1 0 0];
    hand.out{2}.ax.s = 0;
    hand.out{2}.ax.x = fftx;
    hand.out{2}.y = abs(filf_fft(1:sz1/2,:));
    hand.out{2}.title = 'Filtered';
   
end

[fpath,name,ext] = fileparts(get(handles.hh, 'FileName'));
eval([name '(''SetDataSource'', 2, hand)']);

% --------------------------------------------------------------------
function handles = UpdateList(handles)
ListIdx = get(handles.lbParsList, 'Value');

Bands   = get(handles.pmBandSelector, 'String');
BandIdx = get(handles.pmBandSelector, 'Value');

sz = length(handles.Band);
if sz == 0,
  set(handles.lbParsList, 'String',{''},'Value',1);
  set(handles.lbParsList,'Enable', OnOff(0))
  set(handles.slPars,'Enable', OnOff(0))
  set(handles.ePars,'Enable', OnOff(0))

  set(handles.pmBandSelector, 'String','No suppression','Value',1);
  set(handles.pmBandSelector, 'Enable',OnOff(0));
  return;
end

Bands = {};
for ii=1:sz, Bands{end+1} = sprintf('Band: %d',ii); end
if BandIdx > sz, BandIdx = sz; end
set(handles.pmBandSelector, 'String', Bands, 'Value', BandIdx);
set(handles.pmBandSelector, 'Enable',OnOff(1));

List = {};
List{end+1} = sprintf('Pos1: %d',handles.Band(BandIdx).Pos1);
List{end+1} = sprintf('Pos2: %d',handles.Band(BandIdx).Pos2);

if ListIdx > length(List), ListIdx = length(List); end
set(handles.lbParsList, 'String', List ,'Value',ListIdx);
set(handles.lbParsList, 'Enable',OnOff(1));
