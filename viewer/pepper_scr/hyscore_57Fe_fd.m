Sys.S = 0.5
Sys.I = 0.5
Sys.g = [2, 2, 2.06]
Sys.gn = [1; 1]*0.1806
Sys.lw = 0.1
o = '----------------------------'
rotth = 25
Sys.A(1, :) = [-2.4, -2.4, 1.91]
Sys.A(2, :) = [4.68 4.68 2.12]
Sys.Apa(1, :) = [0,30, 30]
Sys.Apa(2, :) = [0,50, 145]
o = '----------------------------'
Exp.gField = 2.
Exp.mwFreq = 34
Exp.ExciteWidth = 50
Exp.tau = 200
Shift.ding = 1
o = '----------------------------'
Opt.Symmetry = 'Ci'
Opt.nKnots = 30
Opt.LineStyle = 'none'
Opt.Sim = 'fd'
Opt.Marker = '.'
Opt.ShowCor = [1, 1; 2,2]
Opt.Verbosity = 2
Opt.KillNeg = 1
Opt.ProdRule = 1
Shift.contour = [0:0.05:1]
Shift.x1 = 0
Shift.x2 = 0
Shift.Scale = 1