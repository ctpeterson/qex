import qex, gauge, physics/qcdTypes, maths/groupOps, gauge/stoutsmear
import math, os, sequtils, strformat, strutils, times, typetraits, json
import purestout, utilityprocs

const 
  nc {.intDefine.} = getDefaultNc()
  dA {.intDefine.} = nc*nc-1

type 
  GaugeActType = enum 
    ActWilson, 
    ActAdjoint, 
    ActRect, 
    ActSymanzik, 
    ActIwasaki, 
    ActDBW2
  NambuActType = enum
    NActWilson, 
    NActAdjoint, 
    NActRect, 
    NActSymanzik, 
    NActIwasaki, 
    NActDBW2,
    NActWilsonStout,  
    NActRectStout, 
    NActSymanzikStout, 
    NActIwasakiStout, 
    NActDBW2Stout,
    NActPlaquetteCharge,
    NActPlaquetteChargeStout

type DAdjointVectorV = Color[VectorArray[dA,DComplexV]]
var Ta = newSeq[Color[MatrixArray[nc,nc,DComplexV]]](dA)

for a in 0..<dA:
  forO b, 0, dA-1:
    forO c, 0, dA-1: Ta[a][b,c] := sugen(nc)[a][b,c]

converter toGaugeActType(s:string): GaugeActType = parseEnum[GaugeActType](s)
converter toNambuActType(s:string): NambuActType = parseEnum[NambuActType](s)

proc AdjointColorVector(l: Layout): auto = l.newField(DAdjointVectorV)

proc newAdjoint(l: Layout): auto =
  result = newSeq[type(l.AdjointColorVector)]()
  for mu in 0..<l.nDim: result.add l.AdjointColorVector

proc projectTa(u: auto; U: auto) =
  threads:
    for mu in 0..<u.len:
      for s in u[mu]:
        forO a, 0, dA-1: u[mu][s][a] := -2.0*trace(Ta[a]*U[mu][s])

proc expandTa(u: auto; v: auto) =
  threads:
    for mu in 0..<v.len:
      for s in v[mu]:
        u[mu][s] := 0
        forO a, 0, dA-1: u[mu][s] += v[mu][s][a]*Ta[a]

qexinit()

letParam:
  gact: GaugeActType = "ActSymanzik"
  nact: NambuActType = "NActPlaquetteCharge"
  lat = @[8,8,8,8]
  beta = 9.0
  nBeta = 9.0
  betaQ = 9.0
  adjFac = -0.25
  nAdjFac = -0.25
  rectFac = -1.0/12.0
  nRectFac = -1.0/12.0
  tau = 1.0
  trajs = 50
  steps = 100
  stepSize = tau/steps
  seed: uint64 = 987654321
  nstout = 3
  rho = 0.1
  
var  
  maxFlowTime = (0.5*lat[^1])*(0.5*lat[^1])/8.0
  flowInfo = %* {
    "C0p0": {
      "action": "Wilson",
      "path": "./",
      "time-increments": [0.02,0.1],
      "maximum-flow-times": [5.0,maxFlowTime]
    }
  }

let
  gc = case gact
    of ActWilson: GaugeActionCoeffs(plaq: beta)
    of ActAdjoint: GaugeActionCoeffs(plaq: beta, adjplaq: beta*adjFac)
    of ActRect: gaugeActRect(beta, rectFac)
    of ActSymanzik: Symanzik(beta)
    of ActIwasaki: Iwasaki(beta)
    of ActDBW2: DBW2(beta)
  ngc = case nact
    of NActWilson,NActWilsonStout: GaugeActionCoeffs(plaq: nBeta)
    of NActAdjoint: GaugeActionCoeffs(plaq: nBeta, adjplaq: nBeta*nAdjFac)
    of NActRect,NActRectStout: gaugeActRect(nBeta, nRectFac)
    of NActSymanzik,NActSymanzikStout: Symanzik(nBeta)
    of NActIwasaki,NActIwasakiStout: Iwasaki(nBeta)
    of NActDBW2,NActDBW2Stout: DBW2(nBeta)
    of NActPlaquetteCharge,NActPlaquetteChargeStout: 
      GaugeActionCoeffs(plaq: betaQ)

var
  lo = lat.newLayout
  vol = lo.physVol

var 
  r = lo.newRNGField(MRG32k3a, seed)
  g = lo.newGauge()
  bg = lo.newGauge()
  f = lo.newGauge()
  bf = lo.newGauge()
  fi = lo.newAdjoint()
  pi = lo.newAdjoint()
  qi = lo.newAdjoint()
  R: MRG32k3a 

var
  stout = case nact
    of NActWilsonStout,NActRectStout: lo.newStoutLinks(nstout,rho)
    of NActSymanzikStout,NActIwasakiStout: lo.newStoutLinks(nstout,rho)
    of NActDBW2Stout: lo.newStoutLinks(nstout,rho)
    else: lo.newStoutLinks(0,0.0)

R.seed(seed, 987654321)

proc action(): float =
  result = case gact
    of ActAdjoint: gc.actionA(g)
    else: gc.gaugeAction1(g)

proc auxAction(): float =
  result = case nact
    of NActAdjoint: ngc.actionA(g)
    of NActPlaquetteCharge: ngc.topologicalCharge(g)
    of NActWilsonStout,NActRectStout: ngc.smearedAction1(stout,g)
    of NActSymanzikStout,NActIwasakiStout: ngc.smearedAction1(stout,g)
    of NActDBW2Stout: ngc.smearedAction1(stout,g)
    of NActPlaquetteChargeStout: 
      stout.smear(g)
      ngc.topologicalCharge(stout.su[^1])
    else: ngc.gaugeAction1(g)

proc updateMomentum(ui: auto; vi: auto; dtau: float) =
  threads:
    for mu in 0..<ui.len:
      for s in ui[mu]:
        forO a, 0, dA-1: 
          ui[mu][s][a] += dtau*fi[mu][s][a]*vi[mu][s][a]

proc nambuForce() =
  case nact:
    of NActAdjoint: ngc.forceA(g,f)
    of NActPlaquetteCharge: ngc.topologicalForce(g,f)
    of NActWilsonStout,NActRectStout: ngc.smearedGaugeForce(stout,g,f)
    of NActSymanzikStout,NActIwasakiStout: ngc.smearedGaugeForce(stout,g,f)
    of NActDBW2Stout: ngc.smearedGaugeForce(stout,g,f)
    of NActPlaquetteChargeStout: 
      stout.smear(g)
      for idx in countdown(stout.nstout,0):
        let xdi = stout.nstout - idx
        if xdi == 0: ngc.topologicalDeriv(stout.su[^1],stout.sf[^1])
        elif xdi == stout.nstout: stout.stout[idx].smearDeriv(f,stout.sf[1])
        else: stout.stout[idx].smearDeriv(stout.sf[idx],stout.sf[idx+1])
      contractProjectTAH(g,f)
      #threads:
      #  for mu in 0..<f.len: f[mu].projectTAH(f[mu])
    else: ngc.gaugeForce(g,f)
  bf.setGauge(f)

proc updateP(dtau: float) =
  fi.projectTa(f)
  pi.updateMomentum(qi,-dtau)

proc updateQ(dtau: float) = 
  fi.projectTa(bf)
  qi.updateMomentum(pi,-dtau)

proc updateU(dtau: float) =
  # Gauge field update
  threads:
    for mu in 0..<fi.len:
      for s in fi[mu]:
         forO a, 0, dA-1: fi[mu][s][a] := pi[mu][s][a]*qi[mu][s][a]
  f.expandTa(fi)
  threads:
    for mu in 0..<g.len:
      for s in g[mu]: g[mu][s] := exp(dtau*f[mu][s])*g[mu][s]

  # Calculate force
  nambuForce()
  case gact:
    of ActAdjoint: gc.forceA(g,f)
    else: gc.gaugeForce(g,f)
  threads:
    for mu in 0..<f.len:
      for s in f[mu]:
        f[mu][s] -= bf[mu][s]

proc drawMomentum(ui: auto) =
  threads: f.randomTAH(r)
  ui.projectTa(f)

proc uiNorm2(ui: auto): float =
  result = 0.0
  for mu in 0..<ui.len: result += ui[mu].norm2
  result *= 0.5

proc hamiltonian(): float =
  let (pi2,qi2,s) = (0.5*pi.uiNorm2(),0.5*qi.uiNorm2(),action())
  result = pi2 + qi2 + s

proc auxHamiltonian(): float =
  let (qi2,s) = (0.5*qi.uiNorm2(),auxAction())
  result = qi2 + s

for traj in 0..<trajs:
  var hi,hf,ahi,ahf,dh,adh: float

  bg.setGauge(g)

  pi.drawMomentum()
  qi.drawMomentum()
  
  (hi,ahi) = (hamiltonian(),auxHamiltonian())

  for step in 0..<steps:
    #qexLog "step: " & $(step)
    if step == 0: updateU(0.5*stepSize)
    else: updateU(stepSize)
    updateP(0.5*stepSize)
    updateQ(stepSize)
    updateP(0.5*stepSize)
    if step == steps-1: updateU(0.5*stepSize)

  (hf,ahf) = (hamiltonian(),auxHamiltonian())

  (dh,adh) = (hf-hi,ahf-ahi)

  if R.uniform <= exp(-dh):
    echo "ACCEPT: dH: ", dH, " dG: ", adh
    g.reunit
  else:
    echo "REJECT: dH: ", dH, " dG: ", adh
    g.setGauge(bg)

  g.mplaq
  g.ploop

  #[
  let filename = "C0p0_" & $(traj) & ".log"
  flowInfo["C0p0"].setFilename(filename)
  g.gradientFlow(flowInfo):
    let output = measurements.formatMeasurements(style="KS_nHYP_FA")
    f.write(output & "\n")
  ]#

qexfinalize()