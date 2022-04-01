function y = tripleprocess(y, bl1, bl2)

sz= size(y);
axx = [1:sz(1)].';
axy = [1:sz(2)].';

% phase in first dimension
for kk=1:sz(2)
  % initial baseline correction (for the phase)
  ss = y(:, kk);
  [a,b] = lin_lsf(axx(bl2, 1), ss(bl2, 1));
  % phase correction   
  th = fminsearch(@m_i_n, 0, [], ss-(b*axx+a));
  ss = ss.*exp(-i*th);
  if abs(max(ss)) < abs(min(ss))
      ss = -ss;
  end
  % baseline correction
  [a,b] = lin_lsf(axx(bl1, 1), ss(bl1, 1));
  y(:, kk) = ss-(b*axx+a);
end

% baseline in second dimension
for kk=1:sz(1)
  ss = y(kk, :).';
  [a,b] = lin_lsf(axy(bl1, 1), ss(bl1, 1));
  y(kk, :) = (ss-(b*axy+a)).';
end

y = -y./max(max(-y));

function [a,b] = lin_lsf(t, Data);

n = size(Data,1);
My = sum(Data)/n;
Mx = sum(t)/n;

b = sum((t-Mx).*(Data-My))/sum((t-Mx).*(t-Mx));
a = My - b * Mx;

function s=m_i_n(th,Data)
s = -sqrt(sum(real(Data*exp(-i*th)).^2));