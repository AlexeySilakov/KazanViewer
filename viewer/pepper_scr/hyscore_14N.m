Sys.S = 0.5
Sys.I = 1
Sys.g(1) = 2
Sys.g(2) = 2.04
Sys.g(3) = 2.1
Sys.gn = 0.4037
Sys.lw = 0.15
o = '----------------------------'
Sys.A(1) = 3.8
Sys.A(2) = 3.5
Sys.A(3) = -0.2
Sys.Apa(1) = 0
Sys.Apa(2) = 22
Sys.Apa(3) = 0
Sys.K = 0.96
Sys.eta = -0.3
Sys.Q(1:3) = [1-Sys.eta, 1+Sys.eta, -2]*Sys.K
Sys.Qpa(1) = 0
Sys.Qpa(2) = 5
Sys.Qpa(3) = 0
o = '----------------------------'
Exp.Field = 324.8
%Exp.mwFreq = 33.865  
Exp.mwFreq = 9.548965
Exp.MaxFreq = 10
Exp.ExciteWidth = 120
Exp.tau = 180
Opt.Symmetry = 'Ci'
Exp.nPoints = 512
Opt.nKnots = 40
Shift.ding = 0
o = '----------------------------'
Opt.LineStyle = '-'
Opt.Sim = 'fd'
Opt.Marker = '.'
Opt.ShowCor = 0
Opt.KillNeg = 1
Shift.contour = [0.0:0.05:1]
Shift.x1 = 0
Shift.x2 = 0
Shift.Scale = 1
Shift.Overlay = 0
Opt.Treshold = 0.2
