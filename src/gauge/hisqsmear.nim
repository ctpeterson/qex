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
    dsdx_dxdw.fat7lDeriv(w,dsdsu,fat7l2,w,dsdsul,naik,info) # Second fat7
    threads: # Unitary projection
      for mu in 0..<dsdx_dxdw_dwdv.len:
        for s in dsdx_dxdw_dwdv[mu]:
          dsdx_dxdw_dwdv[mu][s].projectUderiv(w[mu][s],v[mu][s],dsdx_dxdw[mu][s])
    dsdu.fat7lDeriv(u,dsdx_dxdw_dwdv,fat7l1,info) # First fat7

  if displayPerformance: echo $(info)
  return smearedForce

when isMainModule:
  import gauge/gaugeAction, strformat
  qexInit()
  let
    defaultLat = @[8,8,8,8]
  var
    (lo, g, r) = setupLattice(defaultLat)
    dg = lo.newGauge()
    g2 = lo.newGauge()
    fd = lo.newGauge()
    fl = lo.newGauge()
    ll = lo.newGauge()
    fl2 = lo.newGauge()
    ll2 = lo.newGauge()
    fc = lo.newGauge()
    lc = lo.newGauge()
    ch = lo.newGauge()
    chl = lo.newGauge()
    gc = GaugeActionCoeffs(plaq:1.0)
    eps = floatParam("eps", 1e-6)
    #warm = floatParam("warm", 1e-5)
    hisq = newHISQ()

  #g.unit
  #g.gaussian r
  g.random r
  #g.warm warm, r
  g.stagPhase
  dg.gaussian r
  ch.gaussian r
  chl.gaussian r
  for mu in 0..<g.len:
    for e in g[mu]:
      dg[mu][e] *= eps
      g2[mu][e] := g[mu][e] + dg[mu][e]

  proc check(da: float, tol: float) =
    var ds = 0.0
    for mu in 0..3:
      ds += redot(dg[mu], fd[mu])
    let r = (da-ds)/da
    echo &"  da {da}  ds {ds}  rel {r}"
    if abs(r) > tol*eps:
      echo &"> ERROR rel error |{r}| > {tol*eps}"

  proc checkG =
    #echo "Checking GaugeDeriv"
    #for mu in 0..<fd.len:
    #  fd[mu] := 0
    #let a = gc.gaugeAction2(g)
    #let a2 = gc.gaugeAction2(g2)
    #gc.gaugeDeriv2(g, fd)
    echo "Checking redot deriv"
    var a, a2 = 0.0
    for mu in 0..<fd.len:
      a += redot(ch[mu], g[mu])
      a2 += redot(ch[mu], g2[mu])
      fd[mu] := ch[mu]
    check(a2-a, 1)
  checkG()

  proc checkL(name: string, tol: float) =
    echo "Checking ", name
    resetTimers()
    let f = hisq.smearGetForce(g, fl, ll)
    let f2 = hisq.smearGetForce(g2, fl2, ll2)
    var a, a2 = 0.0
    for mu in 0..<fd.len:
      a += redot(ch[mu], fl[mu]) + redot(chl[mu], ll[mu])
      a2 += redot(ch[mu], fl2[mu]) + redot(chl[mu], ll2[mu])
      fd[mu] := 0
      fc[mu] := ch[mu]
      lc[mu] := chl[mu]
    #gc.gaugeDeriv2(fl, fc)
    #gc.gaugeDeriv2(ll, lc)
    #let a = gc.gaugeAction2(fl) + gc.gaugeAction2(ll)
    #let a2 = gc.gaugeAction2(fl2) + gc.gaugeAction2(ll2)
    f(fd, fc, lc)
    f2(fd, fc, lc)  # combine forces for better accuracy
    check(2.0*(a2-a), tol)  # 2.0 due to combined forces

  hisq.fat7first.oneLink = 1.0
  hisq.fat7first.threeStaple = 0.0
  hisq.fat7first.fiveStaple = 0.0
  hisq.fat7first.sevenStaple = 0.0
  hisq.fat7first.lepage = 0.0
  hisq.fat7second.oneLink = 1.0
  hisq.fat7second.threeStaple = 0.0
  hisq.fat7second.fiveStaple = 0.0
  hisq.fat7second.sevenStaple = 0.0
  hisq.fat7second.lepage = 0.0
  hisq.naik = 0.0
  checkL("oneLink", 1)

  hisq.fat7first.setHisqFat7(0.0, 0.0)
  hisq.fat7second.setHisqFat7(0.0, 0.0)
  hisq.naik = 0.0
  checkL("all", 1)

  hisq = newHisq()
  checkL("all", 1)

  echoProf()
  qexFinalize()
