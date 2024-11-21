import qex
import sequtils, json, strutils, io, sets
import gauge
import gauge/gaugeUtils

proc norm2*(p:auto): float =
  var output: float
  threads:
    var temp = 0.0
    for mu in 0..<p.len: temp += p[mu].norm2
    threadBarrier()
    threadMaster: output = temp
  output

proc reunit*(g:auto) =
  tic()
  threads:
    let d = g.checkSU
    threadBarrier()
    echo "unitary deviation avg: ",d.avg," max: ",d.max
    g.projectSU
    threadBarrier()
    let dd = g.checkSU
    echo "new unitary deviation avg: ",dd.avg," max: ",dd.max
  toc("reunit")

proc mplaq*(g:auto) =
  tic()
  let
    pl = g.plaq
    nl = pl.len div 2
    ps = pl[0..<nl].sum * 2.0
    pt = pl[nl..^1].sum * 2.0
  echo "MEASplaq ss: ",ps,"  st: ",pt,"  tot: ",0.5*(ps+pt)
  toc("plaq")

proc ploop*(g:auto) =
  tic()
  let pg = g[0].l.physGeom
  var pl = newseq[typeof(g.wline @[1])](pg.len)
  for i in 0..<pg.len:
    pl[i] = g.wline repeat(i+1, pg[i])
  let
    pls = pl[0..^2].sum / float(pl.len-1)
    plt = pl[^1]
  echo "MEASploop spatial: ",pls.re," ",pls.im," temporal: ",plt.re," ",plt.im
  toc("ploop")

proc setGauge*(u: auto; v: auto) =
  threads:
    for mu in 0..<u.len:
      for s in u[mu]: u[mu][s] := v[mu][s]

template gradientFlow*(u: auto; info: JsonNode; body: untyped) =
  var 
    v = u[0].l.newGauge
    f {.inject,used.}: File
    tau {.inject.}: float
    measurements {.inject.}: JsonNode
  for flow,flowInfo in info:
    let
      loops = case info[flow].hasKey("loops")
        of true: info[flow]["loops"].getInt()
        of false: 1
      beta = case info[flow].hasKey("beta")
        of true: info[flow]["beta"].getFloat()
        of false: 1.0
      cr = case info[flow]["action"].getStr()
        of "Rectangle": info[flow]["cr"].getFloat()
        else: 0.0
      gc = case info[flow]["action"].getStr()
        of "Rectangle": gaugeActRect(beta,cr)
        else: GaugeActionCoeffs(plaq:beta)
      fn = info[flow]["path"].getStr() & info[flow]["filename"].getStr()
      dts = info[flow].get("time-increments")
      maxFlts = info[flow].get("maximum-flow-times")
    var lastMaxFlt = 0.0
    threads: v := u
    v.reunit
    f = fn.open(fmWrite)
    for (dt,maxFlt) in zip(dts,maxFlts):
      gc.gradientFlow(v,dt):
        tau = flowTime + lastMaxFlt
        measurements = v.flowMeasurements(loops,tau)
        body
        if flowTime>=maxFlt: break
      lastMaxFlt = maxFlt
    f.close()
      
proc formatMeasurements*(
    measurements: JsonNode;
    style: string = "default"
  ): string =
  result = ""
  let tau = measurements["flow-time"].getFloat()
  case style:
    of "default":
      var concatenation = @["FLOW",tau.formatFloat(ffDecimal,3)]
      result = concatenation.join(" ")
      for key,measurement in measurements:
        let O = measurement.getFloat()
        concatenation = @[result,O.formatFloat(ffDecimal,13)] 
        result = concatenation.join(" ")
    of "KS_nHYP_FA": # Style of https://github.com/daschaich/KS_nHYP_FA
      let
        standIn = %* {"stand-in": 0.0}
        observables = @[
          measurements["plaquette"],
          measurements["clover"],
          measurements["t2E-clover"],
          standIn["stand-in"],
          measurements["t2E-plaquette"],
          measurements["topological-charge"],
          measurements["t2E-spacelike-clover"],
          measurements["t2E-timelike-clover"],
          measurements["Re(timelike-Polyakov-loop)"],
          measurements["Im(timelike-Polyakov-loop)"],
          measurements["Re(spacelike-Polyakov-loop)"],
          measurements["Im(spacelike-Polyakov-loop)"]
        ]
      var concatenation = @["FLOW",tau.formatFloat(ffDecimal,2)]
      result = concatenation.join(" ")
      for observable in observables:
        let O = observable.getFloat()
        concatenation = @[result,O.formatFloat(ffDecimal,13)] 
        result = concatenation.join(" ")
    else: discard

proc setFilename*(info: auto; fn: string) =
  info["filename"] = %* fn 