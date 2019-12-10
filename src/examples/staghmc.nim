import qex
import gauge, physics/qcdTypes
import physics/stagSolve
import mdevolve

const ReversibilityCheck {.booldefine.} = false

qexinit()

let
  lat = intSeqParam("lat", @[8,8,8,8])
  #lat = @[8,8,8]
  #lat = @[32,32]
  #lat = @[1024,1024]
  lo = lat.newLayout
  #gc = GaugeActionCoeffs(plaq:6)
  gc = GaugeActionCoeffs(plaq:6,adjplaq:1)
var r = lo.newRNGField(RngMilc6, 987654321)
var R:RngMilc6  # global RNG
R.seed(987654321, 987654321)

var g = lo.newgauge
#g.random r
g.unit

echo 6.0*g.plaq
echo g.gaugeAction2 gc
echo gc.actionA g

var
  p = lo.newgauge
  f = lo.newgauge
  g0 = lo.newgauge
  phi = lo.ColorVector()
  psi = lo.ColorVector()

let mass = floatParam("mass", 0.1)
let stag = newStag(g)
var spa = initSolverParams()
#spa.subsetName = "even"
spa.r2req = floatParam("arsq", 1e-20)
spa.maxits = 10000
var spf = initSolverParams()
#spf.subsetName = "even"
spf.r2req = floatParam("frsq", 1e-12)
spf.maxits = 10000
spf.verbosity = 0

let
  tau = floatParam("tau", 1.0)
  gsteps = intParam("gsteps", 4)
  fsteps = intParam("fsteps", 4)
  trajs = intParam("ntraj", 10)

template rephase(g: typed) =
  g.setBC
  threadBarrier()
  g.stagPhase

proc olf(f: var any, v1: any, v2: any) =
  var t {.noInit.}: type(f)
  for i in 0..<v1.len:
    for j in 0..<v2.len:
      t[i,j] := v1[i] * v2[j].adj
  projectTAH(f, t)

proc oneLinkForce(f: any, p: any, g: any) =
  let t = newTransporters(g, p, 1)
  for mu in 0..<g.len:
    discard t[mu] ^* p
  for mu in 0..<g.len:
    for i in f[mu]:
      olf(f[mu][i], p[i], t[mu].field[i])
    for i in f[mu].odd:
      f[mu][i] *= -1

proc fforce(f: any) =
  tic()
  threads:
    g.rephase
  toc("fforce rephase")
  stag.solve(psi, phi, mass, spf)
  toc("fforce solve")
  #stagD(stag.so, psi, g, psi, 0.0)
  f.oneLinkForce(psi, g)
  toc("fforce olf")
  threads:
    g.rephase
  toc("fforce rephase 2")

proc mdt(t: float) =
  tic()
  threads:
    for mu in 0..<g.len:
      for s in g[mu]:
        g[mu][s] := exp(t*p[mu][s])*g[mu][s]
  toc("mdt")
proc mdv(t: float) =
  tic()
  gc.forceA(g, f)
  threads:
    for mu in 0..<f.len:
      p[mu] -= t*f[mu]
  toc("mdv")

proc mdvf(t: float) =
  tic()
  #let s = t*floatParam("s", 1.0)
  let s = -0.5*t/mass
  f.fforce()
  threads:
    for mu in 0..<f.len:
      p[mu] -= s*f[mu]
  toc("mdvf")

proc mdvf2(t: float) =
  mdv(t)
  mdvf(t)

# For force gradient update
const useFG = true
const useApproxFG2 = false
proc fgv(t: float) =
  tic()
  gc.forceA(g, f)
  threads:
    for mu in 0..<g.len:
      for s in g[mu]:
        g[mu][s] := exp((-t)*f[mu][s])*g[mu][s]
  toc("fgv")
proc fgvf(t: float) =
  tic()
  let t = -0.5*t/mass
  f.fforce()
  threads:
    for mu in 0..<g.len:
      for s in g[mu]:
        g[mu][s] := exp((-t)*f[mu][s])*g[mu][s]
  toc("fgvf")
var gg = lo.newgauge
proc fgsave =
  threads:
    for mu in 0..<g.len:
      gg[mu] := g[mu]
proc fgload =
  threads:
    for mu in 0..<g.len:
      g[mu] := gg[mu]

# Compined update for sharing computations
proc mdvAll(t: openarray[float]) =
  # TODO: actually share computation.
  # For now, just do it separately.
  if t[0] != 0: mdv t[0]
  if t[1] != 0: mdvf t[1]
proc mdvAllfga(ts,gs:openarray[float]) =
  # TODO: actually share computation.
  # For now, just do it separately.
  let
    gt = ts[0] # 0 for gauge
    gg = gs[0]
    ft = ts[1] # 1 for fermion
    fg = gs[1]
  # echo "mdvAll: gauge: ",gt," ",gg,"  fermion: ",ft," ",fg
  # For gauge
  if gg != 0:
    if gt != 0:
      fgsave()
      if useApproxFG2:
        # Approximate the force gradient update with two Taylor expansions.
        let (tf,tg) = approximateFGcoeff2(gt,gg)
        fgv tg[0]
        mdv tf[0]
        fgload()
        fgv tg[1]
        mdv tf[1]
      else:
        # Approximate the force gradient update with a Taylor expansion.
        let (tf,tg) = approximateFGcoeff(gt,gg)
        # echo "gauge fg: ",tf," ",tg
        fgv tg
        mdv tf
      fgload()
    else:
      quit("Force gradient without the force update.")
  elif gt != 0:
    mdv gt
  # For fermion
  if fg != 0:
    if ft != 0:
      fgsave()
      if useApproxFG2:
        # Approximate the force gradient update with two Taylor expansions.
        let (tf,tg) = approximateFGcoeff2(ft,fg)
        fgvf tg[0]
        mdvf tf[0]
        fgload()
        fgvf tg[1]
        mdvf tf[1]
      else:
        # Approximate the force gradient update with a Taylor expansion.
        let (tf,tg) = approximateFGcoeff(ft,fg)
        # echo "fermion fg: ",tf," ",tg
        fgvf tg
        mdvf tf
      fgload()
    else:
      quit("Force gradient without the force update.")
  elif ft != 0:
    mdvf ft

#[ Nested integrators
let
  (V,T) = newIntegratorPair(mdvAll, mdt)
  Hg = mkOmelyan2MN(steps = gsteps div fsteps, V = V[0], T = T)
  H = mkOmelyan2MN(steps = fsteps, V = V[1], T = Hg)
]#
let
  # Omelyan's triple star integrators, see Omelyan et. al. (2003)
  H =
    when useFG:
      let
        (VAllG,T) = newIntegratorPair(mdvAllfga, mdt)
        V = VAllG[0]
        Vf = VAllG[1]
      newParallelEvolution(
        # mkOmelyan4MN4F2GVG(steps = gsteps, V = V, T = T),
        # mkOmelyan4MN4F2GV(steps = gsteps, V = V, T = T),
        # mkOmelyan4MN5F1GV(steps = gsteps, V = V, T = T),
        # mkOmelyan4MN5F1GP(steps = gsteps, V = V, T = T),
        # mkOmelyan4MN5F2GV(steps = gsteps, V = V, T = T),
        mkOmelyan4MN5F2GP(steps = gsteps, V = V, T = T),
        # mkOmelyan6MN5F3GP(steps = gsteps, V = V, T = T),
        mkOmelyan4MN5F2GP(steps = fsteps, V = Vf, T = T))
    else:
      let
        (VAll,T) = newIntegratorPair(mdvAll, mdt)
        V = VAll[0]
        Vf = VAll[1]
      newParallelEvolution(
        # mkOmelyan2MN(steps = gsteps, V = V, T = T),
        # mkOmelyan4MN5FP(steps = gsteps, V = V, T = T),
        # mkOmelyan4MN5FV(steps = gsteps, V = V, T = T),
        mkOmelyan6MN7FV(steps = gsteps, V = V, T = T),
        mkOmelyan6MN7FV(steps = fsteps, V = Vf, T = T))

for n in 1..trajs:
  tic()
  var p2 = 0.0
  var f2 = 0.0
  threads:
    p.randomTAH r
    var p2t = 0.0
    for i in 0..<p.len:
      p2t += p[i].norm2
      g0[i] := g[i]
    threadMaster: p2 = p2t
    psi.gaussian r
    g.rephase
    threadBarrier()
    stag.D(phi, psi, mass)
    threadBarrier()
    phi.odd := 0
  toc("init traj")
  stag.solve(psi, phi, mass, spa)
  toc("fa solve 1")
  threads:
    var psi2 = psi.norm2()
    threadMaster: f2 = psi2
    g.rephase
  let
    ga0 = gc.actionA g0
    fa0 = 0.5*f2
    t0 = 0.5*p2
    h0 = ga0 + fa0 + t0
  toc("init gauge action")
  echo "Begin H: ",h0,"  Sg: ",ga0,"  Sf: ",fa0,"  T: ",t0

  H.evolve tau
  H.finish
  toc("evolve")

  threads:
    var p2t = 0.0
    for i in 0..<p.len:
      p2t += p[i].norm2
    threadMaster: p2 = p2t
    g.rephase
  toc("p norm2, rephase")
  stag.solve(psi, phi, mass, spa)
  toc("fa solve 2")
  threads:
    var psi2 = psi.norm2()
    threadMaster: f2 = psi2
    g.rephase
  let
    ga1 = gc.actionA g
    fa1 = 0.5*f2
    t1 = 0.5*p2
    h1 = ga1 + fa1 + t1
  toc("final gauge action")
  echo "End H: ",h1,"  Sg: ",ga1,"  Sf: ",fa1,"  T: ",t1

  when ReversibilityCheck:
    block:
      var g1 = lo.newgauge
      var p1 = lo.newgauge
      threads:
        for i in 0..<g1.len:
          g1[i] := g[i]
          p1[i] := p[i]
          p[i] := -1*p[i]
      H.evolve tau
      H.finish
      threads:
        var p2t = 0.0
        for i in 0..<p.len:
          p2t += p[i].norm2
        threadMaster: p2 = p2t
        g.rephase
      toc("p norm2, rephase")
      stag.solve(psi, phi, mass, spa)
      toc("fa solve 2")
      threads:
        var psi2 = psi.norm2()
        threadMaster: f2 = psi2
        g.rephase
      let
        ga1 = gc.actionA g
        fa1 = 0.5*f2
        t1 = 0.5*p2
        h1 = ga1 + fa1 + t1
      echo "Reversed H: ",h1,"  Sg: ",ga1,"  Sf: ",fa1,"  T: ",t1
      echo "Reversibility: dH: ",h1-h0,"  dSg: ",ga1-ga0,"  dSf: ",fa1-fa0,"  dT: ",t1-t0
      #echo p[0][0]
      for i in 0..<g1.len:
        g[i] := g1[i]
        p[i] := p1[i]
    toc("reversibility")

  let
    dH = h1 - h0
    acc = exp(-dH)
    accr = R.uniform
  if accr <= acc:  # accept
    echo "ACCEPT:  dH: ",dH,"  exp(-dH): ",acc,"  r: ",accr
  else:  # reject
    echo "REJECT:  dH: ",dH,"  exp(-dH): ",acc,"  r: ",accr
    threads:
      for i in 0..<g.len:
        g[i] := g0[i]

  echo 6.0*g.plaq

echoTimers()
qexfinalize()
