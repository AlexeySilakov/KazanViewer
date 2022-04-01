% out = kv_expdata(t1, t2, xo, x, y)
%   Returns  out =  A*exp(-1*(x-xo)/t1) + B*exp(-1*(x-xo)/t2) + C
%   A, B, C are calculated according to [x, y] by LSQ
% alsi 25-Mar-2004

function out = kv_expdata(t1, t2, xo, x, y)
[ss, ind] = max(size(x));
if ind == 2, x = x.'; end
[ss, ind] = max(size(y));
if ind == 2, y = y.'; end

F1 = exp(-1*(x-xo)/t1);
F2 = exp(-1*(x-xo)/t2);

A1 = F1.'*F1; % sum(F1.^2)
B1 = F2.'*F1; %sum(F2.*F1);
C1 = sum(F1);
D1 = F1.'*y;

A2 = F1.'*F2; % sum(F1.^2)      {A1 B1 C1}   {D1}
B2 = F2.'*F2; %sum(F2.*F1);     {A2 B2 C2} = {D2} 
C2 = sum(F2); %                 {A3 B3 C3}   {D3}
D2 = F2.'*y;

A3 = sum(F1); 
B3 = sum(F2); %sum(F2.*F1);
C3 = size(x, 1); %sum(1)
D3 = sum(y);

D = det([A1, B1, C1; A2, B2, C2; A3, B3, C3]);
A = det([D1, B1, C1; D2, B2, C2; D3, B3, C3])/D;
B = det([A1, D1, C1; A2, D2, C2; A3, D3, C3])/D;
C = det([A1, B1, D1; A2, B2, D2; A3, B3, D3])/D;
     
out = A*F1 + B*F2 + C;
if ind == 2, out = out.'; end