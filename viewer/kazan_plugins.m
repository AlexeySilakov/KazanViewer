function out = kazan_plugins
% This is an internal function for kazan viewer
% that tells it where to look for plugins
[fpath,n,e] = fileparts(which('kazan'));
plugins_path = [fpath , filesep, 'plugins', filesep];

out = {};
out{end+1}.name = 'File Description';
out{end}.file = 'description'; % WARNING: it is case sensetive

out{end+1}.name = 'Compare data(compbox)';
out{end}.file = 'CompBox'; % WARNING: it is case sensetive

out{end+1}.name = 'Simulate 1D data(simplugin)';
out{end}.file = 'simplugin'; % WARNING: it is case sensetive

out{end+1}.name = 'Fourier Transform';
out{end}.file = 'fftplugin'; % WARNING: it is case sensetive

out{end+1}.name = 'Fitting box';
out{end}.file = 'fittingbox'; % WARNING: it is case sensetive

out{end+1}.name = 'Cut and Math';
out{end}.file = 'mathplugin'; % WARNING: it is case sensetive

out{end+1}.name = 'Data subtraction (minusbox)';
out{end}.file = 'minusbox'; % WARNING: it is case sensetive

out{end+1}.name = 'Simulate 2D data (hyscorebox)';
out{end}.file = 'hyscorebox'; % WARNING: it is case sensetive

out{end+1}.name = 'Manipulate IR data (FTIR box)';
out{end}.file = 'FTIRBox'; % WARNING: it is case sensetive

out{end+1}.name = 'Singular Value Decomposition (svdplugin)';
out{end}.file = 'svdplugin'; % WARNING: it is case sensetive

out{end+1}.name = 'Import Data from Image (image2data)';
out{end}.file = 'image2data'; % WARNING: it is case sensetive

out{end+1}.name = 'DEER Analysis (deerplugin)';
out{end}.file = 'deerplugin'; % WARNING: it is case sensetive