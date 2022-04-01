function varargout = mqmasplugin(varargin)
% MQMASPLUGIN Application M-file for mqmasplugin.fig
%    FIG = MQMASPLUGIN launch mqmasplugin GUI.
%    MQMASPLUGIN('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 04-Dec-2004 11:48:06

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
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
function FileUpdate(handles)


% --------------------------------------------------------------------
function PlotThis(ha)

h = guidata(ha.figure1);
hand = guidata(h.hh);
hand.out = {};

stage = get(ha.pmStage,'Value');

% determine what is a source file
[path,name,ext] = fileparts(hand.src.filename);
if strfind(name,'mas_ae')
    aecho = hand.src.y;
    pos = findstr(hand.lscript{1},'_ae');
    fname = hand.lscript{1}([1:pos,pos+2:end]);
    eval(fname)
    echo  = y; 
elseif strfind(name,'mas_e')
    echo = hand.src.y;
    pos = findstr(hand.lscript{1},'_e');
    fname = hand.lscript{1}([1:pos,pos-2,pos+1:end]);
    eval(fname)
    aecho  = y; 
else return;
end

[nt1,nt2] = size(echo);

if get(ha.cbSwap, 'Value')
   tmp = aecho; aecho = echo; echo = tmp; clear tmp
end

aechom(:,1:nt2) = aecho(:,nt2:-1:1);
dat = aechom;
% append echo to antiecho matrix
dat(:,nt2+1:2*nt2) = echo(:, :);

t1 = hand.src.ax.x; dt1 = (t1(end)-t1(1))/(length(t1)-1);
t2 = hand.src.ax.y; dt2 = (t2(end)-t2(1))/(length(t2)-1);

%  Fourie Transform FT 2
if stage > 1
    clear awin
    width = str2num(get(ha.eApo2,'String'));
    npoints = floor(width/dt2+.5);
    aw = apowin('gau', 2*npoints, 0.6);
    awin = aw(end)*ones(nt2*2,1);
    sz = 2*min(npoints,nt2);
    sh = nt2 - min(sz,nt2);
    awin([1:sz]+sh)=aw(1:sz);
    
    ft2 = [-nt2:nt2-1]'/dt2/(2*nt2-1);
    fdat = fftshift(fft(fftshift(dat.*awin(:,ones(nt1,1))',2), [], 2),2); 
end

%  Shearing
if stage > 2
    clear fshdat
    if get(ha.pmSpin, 'Value') == 1
        sh = -9/7*pi/nt2; % 3/2
    else
        sh = -12/19*pi/nt2; % 5/2
    end
    
    for jj=1:nt1
        fshdat(jj,:) = fdat(jj,:).*exp(-i*sh*(jj-1)*([1:2*nt2]-(nt2+1)));
    end
end

%  Fourie Transform FT 1
if stage > 3
    clear awin
    width = str2num(get(ha.eApo1,'String'));
    npoints = floor(width/dt1+.5);
    aw = apowin('gau+', npoints, 0.6);
    awin = aw(end)*ones(nt1,1);
    awin(1:min(npoints,nt1))=aw(1:min(npoints,nt1));
    
    ft1 = [-nt1:nt1-1]'/dt1/(2*nt1-1);
    ffdat = fftshift(fft(fshdat.*awin(:,ones(2*nt2,1)), nt1*2, 1),1); 
end

switch stage
case 1,
    hand.out{1}.y = dat';
    hand.out{1}.dsc = dsc;
    hand.out{1}.ax=ax;
    hand.out{1}.ax.x = [-t2(end:-1:1);t2];
    hand.out{1}.ax.y = t1;
    hand.out{1}.ax.xlabel =  't2, ms';
    hand.out{1}.ax.ylabel =  't1, ms';
case 2,
    hand.out{1}.y = fdat';
    hand.out{1}.dsc = dsc;
    hand.out{1}.ax=ax;
    hand.out{1}.ax.x = ft2;
    hand.out{1}.ax.y = t1;
    hand.out{1}.ax.xlabel = 't2, kHz';
    hand.out{1}.ax.ylabel = 't1, ms';    
case 3,
    hand.out{1}.y = fshdat';
    hand.out{1}.dsc = dsc;
    hand.out{1}.ax=ax;
    hand.out{1}.ax.x = ft2;
    hand.out{1}.ax.y = t1;
    hand.out{1}.ax.xlabel = 't2, kHz';
    hand.out{1}.ax.ylabel = 't1, ms'; 
case 4,
    hand.out{1}.y = ffdat;
    hand.out{1}.dsc = dsc;
    hand.out{1}.ax=ax;
    hand.out{1}.ax.x = ft1;
    hand.out{1}.ax.y = ft2;
    hand.out{1}.ax.xlabel = 't1, kHz';
    hand.out{1}.ax.ylabel = 't2, kHz';
end

[path,name,ext] = fileparts(get(h.hh, 'FileName'));
eval([name '(''SetDataSource'', 2, hand)']);

% --------------------------------------------------------------------
function varargout = eApo1_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = eApo2_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = pmSpin_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = pbFFT_Callback(h, eventdata, handles, varargin)
PlotThis(handles);

% --------------------------------------------------------------------
function varargout = cbSwap_Callback(h, eventdata, handles, varargin)
PlotThis(handles);

% --------------------------------------------------------------------
function varargout = pm2_Callback(h, eventdata, handles, varargin)
PlotThis(handles);

