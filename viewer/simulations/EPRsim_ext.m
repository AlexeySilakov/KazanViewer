function varargout = EPRsim_ext(varargin)

S = 1/2;
I = 1/2;
g = [2, 2, 2.065];
A = [0, 0, 0];
Apa = -1;
lw =1;
ee = [0, 0, 0];
D = [0, 0, 0];

HStrain = [0, 0, 0];

Range = [300, 400];
mwFreq = 9.4;
nPoints = 1000;
Temperature = 4;
DetAng = 90;
Harmonic = 1;
nKnots = 10;
nExtraKnots= 10;
nRanges = 5;
PrintSystem = 1;
EigenField = 0;
[fpath, ~, ~] = fileparts(which('EPRsim_ext'));

% % kind of global safeget
if nargin>=1
    for kk = 1:nargin
        if isstruct(varargin{kk})
            Datas = struct2cell(varargin{kk});
            Names = fieldnames(varargin{kk});
            for ii = 1:length(Datas)
                if isnumeric(Datas{ii})
                    tmp = Datas{ii};
                    eval([Names{ii}, '= tmp;']);
                else
                    eval([Names{ii}, '= ''',Datas{ii}, ''';']);
                end    
            end
        end
    end
end

fid = fopen([fpath, filesep, 'EPRsim.inp'], 'w');
nSpins = length(S);
nNucs = length(I)*(sum(A(:))~=0)*(sum(I(:))~=0);
fprintf(fid, 'nSpins %d \n', length(S));
fprintf(fid, 'S_values ');
    for ii = 1:length(S)
        fprintf(fid, '%3.2f ', S(ii));
    end
    fprintf(fid, '\n');
    
for ii = 1:length(S)
    fprintf(fid, 'g_values %f %f %f \n', g(ii, 1), g(ii, 2), g(ii, 3));
end

fprintf(fid, 'fwhh %f \n', lw(1)); % only gaussian goes into the program, lorentzian is convoluted afterwards
fprintf(fid, 'HStrain %f %f %f \n', HStrain(1),HStrain(2), HStrain(3));
fprintf(fid, 'nNucs %d \n', nNucs) ;
fprintf(fid, 'I_values ');
    for ii = 1:length(I)
        fprintf(fid, '%3.2f ', I(ii));
    end
    fprintf(fid, '\n');


% conversion from "EasySpin" notation to EPRsim
idx = [1, 4, 5; 6, 2, 7; 8, 9, 3];
tIdx = [];
A_values = zeros(9, 3);
Aang_values = zeros(9, 3);
if size(Apa, 2)==1
    Apa = A*0;
end
if size(Apa)~=size(A), 
    fclose(fid); 
    error('Size of A does not match size of Apa'); 
end
for ii = 1:nSpins
    if (size(A, 1)<ii), break; end
    for jj = 1:nNucs
        if (size(A, 2)<(3*jj)), break; end
        A_values(idx(ii, jj), :) = A(ii, 3*(jj-1)+(1:3));
        Aang_values(idx(ii, jj), :) = Apa(ii, 3*(jj-1)+(1:3));
        tIdx(end+1) = idx(ii, jj);
    end
end
for ii = 1:max(tIdx)
    fprintf(fid, 'A_values %f %f %f \n', A_values(ii, 1), A_values(ii, 2), A_values(ii, 3));
    fprintf(fid, 'Aang_values %f %f %f \n', Aang_values(ii, 1), Aang_values(ii, 2), Aang_values(ii, 3));
end

if nSpins>1
    count = 1;
    for ii = 1:(nSpins-1)
        for jj = (ii+1):nSpins
            fprintf(fid, 'Ee_values %f %f %f \n',  ee(count, 1), ee(count, 2), ee(count, 3));
            count=count+1;
        end
    end
end

if size(D, 2)==1
    D_values = [-1/3*D(:), -1/3*D(:), 2/3*D(:)];
elseif size(D, 2)==2
    D_values = [-1/3*D(:, 1), -1/3*D(:, 1), 2/3*D(:, 1)] + [1*D(:, 2), -1*D(:, 2), 0*D(:, 2)];
else
    D_values = D;
end 
for ii = 1:size(D, 1)
    fprintf(fid, 'D_values %f %f %f \n', D_values(ii, 1), D_values(ii, 2), D_values(ii, 3));
end

% #
% ### ZeroField S1, S2, S3. One angle set per interaction. 
% ### In principal values Dxx, Dyy, Dzz. 
% D_values 0.0 0.0 0.0 
% D_values 0.0 0.0 0.0 
% #
% Dang_values 0.0 0.0 0.0 
% Dang_values 0.0 0.0 0.0 
% #
% ###  EXPERIMENAL AND OPTIONAL PARAMETERS ####
% ### Field values are in [mT] 

fprintf(fid, 'Bmin %f \n', Range(1)); % mT
fprintf(fid, 'Bmax %f \n', Range(2));
fprintf(fid, 'mwFreq %f \n', mwFreq*1e3); %% -microwave frequency in [GHz] -> [MHz] 
fprintf(fid, 'nPoints %i \n', nPoints); 
fprintf(fid, 'Temperature %i \n', Temperature); %[K]
fprintf(fid, 'DetAng %f \n', DetAng); % in degree sections. 
fprintf(fid, 'Harmonic %i \n', Harmonic);
fprintf(fid, 'nKnots %i \n', nKnots);
fprintf(fid, 'nExtraKnots %i \n', nExtraKnots);
fprintf(fid, 'nRanges %i \n', nRanges);
fprintf(fid, 'PrintSystem %i \n', PrintSystem);
fprintf(fid, 'EigenField %i \n', EigenField);
fclose(fid);

aaa = cd;
cd(fpath);
% tic
system(['del out.dat']);
system('EPRsim.exe -a EPRsim.inp > EPRsim.log');
% toc
cd(aaa);

fid = fopen([fpath, filesep, 'out.dat'], 'r');
count = 0;
while fid<0 && count<10
%     pause(0.0001);
    fid = fopen([fpath, filesep, 'out.dat'], 'r');
    count = count+1;
end

if fid<0, error('can not read the output file'); end
try
    A = textscan(fid, '%f', nPoints*2);
end
fclose(fid);

x  =A{1}(1:2:end);
y  =A{1}(2:2:end);


if length(lw)==2
    if lw(2)~=0
        iy = ifft(y);
        lor = lw(2)/2/pi./((x-mean(x)).^2+(lw(2)/2)^2);
        iw = ifft(lor);
        y = real(fftshift(fft(iy.*iw)));
    end
end
if ~nargout
    figure(12); clf; plot(x, y);
else
    varargout{1} = x*10;% mT>G
    varargout{2} = y;
end




