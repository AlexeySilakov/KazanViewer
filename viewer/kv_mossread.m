function varargout = kv_mossread(filename, varargin)
% read ASCII data of the ".moss" files


fid = fopen(filename, 'r');

if fid == -1, error('sorry, file not found'); return; end

%%%% first string -> start point, ???, folding point, npoints
a = fgets(fid);
F = sscanf(a, ['%f %*s %f %f %f']);
units = sscanf(a, ['%*f %s %*f %*f %*f']);
    startpoint = F(1);
    zervel = (F(2));
    foldingpoint  = (F(3));
    npoints  = F(4);
    
    dX = abs(startpoint/(zervel));
%     dX = abs(startpoint/(zervel));
%%%% second string -> description
description = fgets(fid);
%%%%%%%%%%% then data
tA = textscan (fid,'%f', 'MultipleDelimsAsOne', 1);

fclose(fid);

    data = tA{1};
% data = textread(filename, 'headerlines', 2);

ax = [];
y = [];
dsc = [];

ax.xlabel = 'Velosity, mm/s';
ax.type = 'data';
ax.title = description(1:end-2);
fold = (foldingpoint>0);
if ~fold
    y = data;
    ax.x = ([1:(npoints*2)] - 1)';
    dsc.title = description;
else
    ind = round(foldingpoint);
%     y = zeros(npoints, 1);
    len = size(data, 1);
    
    if ind>round(len/2)
        datalen = len-ind-1;
        st1 = ind-datalen; 
        en2 = len;
    else
        datalen = ind-1;
        st1 = 2; 
        en2 = ind*2-1;
    end
    y1 = data(2:ind);
%     x1 = -((2:ind) - zervel-1)*dX;
    x1 = -((2:ind) -foldingpoint + round(zervel)-1 -0.5)*dX;
    
    y2 = data((ind):len);
%     x2 = (((ind):len)-2*foldingpoint + zervel-1)*dX;
    x2 = (((ind):len)-foldingpoint - round(zervel)-1+0.5)*dX;

    st = max([min(x1), min(x2)]);
    en = min([max(x1), max(x2)]);
    axiss = linspace(st, en, 256);
%       figure(3211); plot(x1, y1); hold on; plot(x2, y2);
    y = spline(x1, y1, axiss) + spline(x2, y2, axiss);
    
%     y = data((ind):-1:st1) + data((ind+1):en2);
%     y = spline(1:ind*dX, data(1:ind), linspace(st1, ind, npoints) )
%     y = spline(st1:ind, data((ind):-1:st1), linspace(st1, ind, npoints) ) + spline( (ind+1):en2, data( (ind+1):en2), linspace(ind+1, en2, npoints) );
    
    ax.x = axiss(:); %([1:length(y)]-1)'*dX + startpoint+0.102;
    dsc.title =  description(1:end-2);
end
dsc.foudlingpoint = foldingpoint;
dsc.startpoint = startpoint;
dsc.zervel = zervel;
% dsc.npoints = npoints;
dsc.stem = dX;
varargout{1} = ax;
varargout{2} = y(:);
varargout{3} = dsc;

