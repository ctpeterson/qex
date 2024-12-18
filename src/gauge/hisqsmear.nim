import qex
import physics/[hisqLinks]
import gauge
import gauge/[fat7l,fat7lderiv]

export hisqLinks

proc newHISQ*(lepage: float = 0.0; naik: float = 1.0): HisqCoefs =
  result = HisqCoefs(naik: -naik/24.0)
  result.fat7first.setHisqFat7(lepage,0.0)
  result.fat7second.setHisqFat7(2.0-lepage,naik)

proc smearGetForce*[T](
    self: HisqCoefs;
    u: T;
    su,sul: T;
    displayPerformance: bool = false
  ): proc(dsdu: var T; dsdsu,dsdsul: T) =
  mixin projectU,projectUderiv
  let
    lo = u[0].l
    fat7l1 = self.fat7first
    fat7l2 = self.fat7second
    naik = self.naik
  var
    v = newOneOf(u)
    w = newOneOf(u)
    info: PerfInfo

  # Smear
  v.makeImpLinks(u,fat7l1,info) # First fat7
  threads: # Unitary projection
    for mu in 0..<w.len:
      for s in w[mu]: w[mu][s].projectU(v[mu][s])
  makeImpLinks(su,w,fat7l2,sul,w,naik,info) # Second fat7

  # Chain rule - retains a reference to u,su,sul
  proc smearedForce(dsdu: var T; dsdsu,dsdsul: T) =
    var
      dsdx_dxdw = newOneOf(dsdu)
      dsdx_dxdw_dwdv = newOneOf(dsdu)
    dsdx_dxdw.fat7lDeriv(su,dsdsu,fat7l2,sul,dsdsul,naik,info) # Second fat7
    threads: # Unitary projection
      for mu in 0..<dsdx_dxdw_dwdv.len:
        for s in dsdx_dxdw_dwdv[mu]:
          dsdx_dxdw_dwdv[mu][s].projectUderiv(w[mu][s],v[mu][s],dsdx_dxdw[mu][s])
    dsdu.fat7lDeriv(u,dsdx_dxdw_dwdv,fat7l1,info) # First fat7

  if displayPerformance: echo $(info)
  return smearedForce

if isMainModule:
  qexInit()
  let
    defaultLat = @[8,8,8,8]
    hisq = newHISQ()
  var
    (lo, g, r) = setupLattice(defaultLat)
    sg = lo.newGauge()
    sgl = lo.newGauge()
    f = lo.newGauge()
    ff = lo.newGauge()
    ffl = lo.newGauge()
    g2 = lo.newGauge()
    sg2 = lo.newGauge()
    sgl2 = lo.newGauge()
    dg = lo.newGauge()
    fd = lo.newGauge()
  g.random
  for mu in 0..<dg.len:
    dg[mu] := 0.00001 * g[mu]
    g2[mu] := g[mu] + dg[mu]
    fd[mu] := 0

  echo g.plaq
  echo g2.plaq
  #hisq.smear(g,sg,sgl)
  var force = hisq.smearGetForce(g,sg,sgl)
  f.force(dg,dg)
  hisq.smear(g2,sg2,sgl2)

  echo f.plaq
  echo sg.plaq
  echo sgl.plaq
  echo "--"
  echo sg.plaq
  echo sgl.plaq
  ff.gaugeForce(sg)
  ffl.gaugeForce(sgl)
  echo "--"
  echo f.plaq
  echo sg.plaq
  echo sgl.plaq
  qexFinalize()
