import os
import basicOps
import strUtils
import qex
import qcdTypes
import stdUtils
import times
import cg
import types
import profile
import metaUtils
import gaugeUtils

#{.emit:"#define memset(a,b,c)".}

type StaggeredD*[T] = object
  sf*:seq[ShiftB[T]]
  sb*:seq[ShiftB[T]]
  sub*:string
  subset*:Subset
type Staggered*[G,T] = object
  se*,so*:StaggeredD[T]
  g*:seq[G]

template initStagDT*(l:var Layout; T:typedesc; ss:string):untyped =
  var sd:StaggeredD[T]
  sd.sf.newSeq(4)
  sd.sb.newSeq(4)
  for mu in 0..<4:
    initShiftB(sd.sf[mu], l, T, mu, 1, ss)
    initShiftB(sd.sb[mu], l, T, mu,-1, ss)
  sd.sub = ss
  sd.subset.layoutSubset(l, ss)
  sd

proc initStagD*(x:Field; sub:string):auto =
  result = initStagDT(x.l, type(x[0]), sub)

template initStagD3T*(l:var Layout; T:typedesc; ss:string):untyped =
  var sd:StaggeredD[T]
  sd.sf.newSeq(8)
  sd.sb.newSeq(8)
  for mu in 0..<4:
    initShiftB(sd.sf[2*mu  ], l, T, mu, 1, ss)
    initShiftB(sd.sf[2*mu+1], l, T, mu, 3, ss)
    initShiftB(sd.sb[2*mu  ], l, T, mu,-1, ss)
    initShiftB(sd.sb[2*mu+1], l, T, mu,-3, ss)
  sd.sub = ss
  sd.subset.layoutSubset(l, ss)
  sd

proc initStagD3*(x:Field; sub:string):auto =
  result = initStagD3T(x.l, type(x[0]), sub)

template stagDPN*(sd:openArray[StaggeredD]; r:openArray[Field];
                  g:openArray[Field2]; x:openArray[Field3];
                  expFlops:int; exp:untyped) {.dirty.} =
  tic()
  let n = sd.len
  #[
  for mu in 0..<g.len:
    for i in 0..<n:
      startSB(sd[i].sf[mu], x[i][ix])
  toc("startShiftF")
  for mu in 0..<g.len:
    for i in 0..<n:
      startSB(sd[i].sb[mu], g[mu][ix].adj*x[i][ix])
  toc("startShiftB")
  ]#
  #var rir = newAlignedMem[type(load1(r[0][0]))](n)
  #var rir:array[10,type(load1(r[0][0]))]
  #for ir{.inject.} in r[0][sd[0].subset]:
  #for ir in r[0][sd[0].subset]:
  #  for i in 0..<n:
  let ns0 = sd[0].subset.lowOuter
  let ns1 = sd[0].subset.highOuter
  let ns = ns1 - ns0
  #tFor ir, ns0, <ns1:
  #  for i in 0..<n:
  let nsn = ns*n
  let nr = 8
  let ni = 16
  let nin = ni * n
  let ninr = ni * nr
  let n4 = n div nr
  let n4r = n mod nr
  tFor iri, 0, <nsn:
      let lr = iri mod nin
      let lr64 = lr div ninr
      let lrr = lr - ninr*lr64
      var nn = nr
      if lr64>=n4: nn = n4r
      let ir = ns0 + (iri div nin)*ni + (lrr div nn)
      let i = lr64*nr + (lrr mod nn)
      #let ir = ns0 + (iri div n)
      #let i = iri mod n
      var rir{.inject,noInit.}:type(load1(r[i][ir]))
      exp
      for mu in 0..<g.len:
        #if mu<g.len-1:
        #  prefetchSB(sd[i].sf[mu+1], ir, x[i][ix])
        #else:
        #  prefetchSB(sd[i].sb[0], ir, x[i][ix])
        localSB(sd[i].sf[mu], ir, imadd(rir, g[mu][ir], it),load1(x[i][ix]))
      for mu in 0..<g.len:
        #if mu<g.len-1:
        #  prefetchSB(sd[i].sb[mu+1], ir, x[i][ix])
        #else:
        #  prefetchSB(sd[i].sf[mu], ir, x[i][ix])
        localSB(sd[i].sb[mu], ir, isub(rir, it), g[mu][ix].adj*x[i][ix])
      assign(r[i][ir], rir)
  toc("local", flops=n*(expFlops+g.len*(72+66+6))*sd[0].subset.len)
  #[
  for mu in 0..<g.len:
    for i in 0..<n:
      boundarySB(sd[i].sf[mu], imadd(r[i][ir], g[mu][ir], it))
  toc("boundaryF")
  for mu in 0..<g.len:
    for i in 0..<n:
      boundarySB(sd[i].sb[mu], isub(r[i][ir], it))
  toc("boundaryB")
  ]#
template stagDMN*(sd:openArray[StaggeredD]; r:openArray[Field];
                  g:openArray[Field2]; x:openArray[Field3];
                  expFlops:int; exp:untyped) =
  tic()
  let n = sd.len
  for mu in 0..<g.len:
    for i in 0..<n:
      startSB(sd[i].sf[mu], x[i][ix])
  toc("startShiftF")
  for mu in 0..<g.len:
    for i in 0..<n:
      startSB(sd[i].sb[mu], g[mu][ix].adj*x[i][ix])
  toc("startShiftB")
  var rir{.inject,noInit.}:seq[type(load1(r[0][0]))]
  rir.newSeq(n)
  for ir{.inject.} in r[sd.subset]:
    exp
    for mu in 0..<g.len:
      for i in 0..<n:
        localSB(sd[i].sf[mu], ir, imsub(rir[i], g[mu][ir], it),load1(x[i][ix]))
      for i in 0..<n:
        localSB(sd[i].sb[mu], ir, iadd(rir[i], it), g[mu][ix].adj*x[i][ix])
    for i in 0..<n:
      assign(r[i][ir], rir[i])
  toc("local", flops=(expFlops+g.len*(72+66+6))*sd.subset.len)
  for mu in 0..<g.len:
    for i in 0..<n:
      boundarySB(sd[i].sf[mu], imsub(r[i][ir], g[mu][ir], it))
  toc("boundaryF")
  for mu in 0..<g.len:
    for i in 0..<n:
      boundarySB(sd[i].sb[mu], iadd(r[i][ir], it))
  toc("boundaryB")

template stagDP*(sd:StaggeredD; r:Field; g:openArray[Field2];
                 x:Field3; expFlops:int; exp:untyped) =
  tic()
  for mu in 0..<g.len:
    startSB(sd.sf[mu], x[ix])
  toc("startShiftF")
  for mu in 0..<g.len:
    startSB(sd.sb[mu], g[mu][ix].adj*x[ix])
  toc("startShiftB")
  for ir{.inject.} in r[sd.subset]:
    var rir{.inject,noInit.}:type(load1(r[ir]))
    exp
    for mu in 0..<g.len:
      localSB(sd.sf[mu], ir, imadd(rir, g[mu][ir], it), load1(x[ix]))
    for mu in 0..<g.len:
      localSB(sd.sb[mu], ir, isub(rir, it), g[mu][ix].adj*x[ix])
    assign(r[ir], rir)
  toc("local", flops=(expFlops+g.len*(72+66+6))*sd.subset.len)
  for mu in 0..<g.len:
    boundarySB(sd.sf[mu], imadd(r[ir], g[mu][ir], it))
  toc("boundaryF")
  for mu in 0..<g.len:
    boundarySB(sd.sb[mu], isub(r[ir], it))
  #threadBarrier()
  toc("boundaryB")

template nVecs(x:untyped):untyped =
  when compiles(nrows(x)): nrows(x)
  else: 1
template getVec(x,i:untyped):untyped = row(x,i)
template setVec(r,x,i:untyped):untyped = setRow(r,x,i)
template stagDP2*(sd:StaggeredD; r:Field; g:openArray[Field2];
                  x:Field3; expFlops:int; exp:untyped) =
  tic()
  #[
  for mu in 0..<len(g):
    startSB(sd.sf[mu], x[ix])
  toc("startShiftF")
  for mu in 0..<g.len:
    startSB(sd.sb[mu], g[mu][ix].adj*x[ix])
  toc("startShiftB")
  ]#
  const n = nVecs(x[0])
  echoImm:
    $n & ": " & $len(x[0])
  for ir{.inject.} in r[sd.subset]:
    for ic{.inject.} in 0..<n:
    #forStatic ic, 0, <n:
    #  block:
        var rir{.inject,noInit.}:type(getVec(r[ir],0))
        for mu in 0..<g.len:
          localSB(sd.sf[mu], ir, imadd(rir, g[mu][ir], it), getVec(x[ix],ic))
        for mu in 0..<g.len:
          localSB(sd.sb[mu], ir, isub(rir, it), g[mu][ix].adj*getVec(x[ix],ic))
        setVec(r[ir], rir, ic)
  toc("local", flops=n*(expFlops+g.len*(72+66+6))*sd.subset.len)
  #[
  for mu in 0..<g.len:
    boundarySB(sd.sf[mu], imadd(r[ir], g[mu][ir], it))
  toc("boundaryF")
  for mu in 0..<g.len:
    boundarySB(sd.sb[mu], isub(r[ir], it))
  #threadBarrier()
  toc("boundaryB")
  ]#
template stagDM*(sd:StaggeredD; r:Field; g:openArray[Field2];
                 x:Field3; expFlops:int; exp:untyped) =
  tic()
  for mu in 0..<g.len:
    startSB(sd.sf[mu], x[ix])
  toc("startShiftF")
  for mu in 0..<g.len:
    startSB(sd.sb[mu], g[mu][ix].adj*x[ix])
  toc("startShiftB")
  for ir{.inject.} in r[sd.subset]:
    var rir{.inject,noInit.}:type(load1(r[ir]))
    exp
    for mu in 0..<g.len:
      localSB(sd.sf[mu], ir, imsub(rir, g[mu][ir], it), load1(x[ix]))
      localSB(sd.sb[mu], ir, iadd(rir, it), g[mu][ix].adj*x[ix])
    assign(r[ir], rir)
  toc("local", flops=(expFlops+g.len*(72+66+6))*sd.subset.len)
  for mu in 0..<g.len:
    boundarySB(sd.sf[mu], imsub(r[ir], g[mu][ir], it))
  toc("boundaryF")
  for mu in 0..<g.len:
    boundarySB(sd.sb[mu], iadd(r[ir], it))
  #threadBarrier()
  toc("boundaryB")

# r = a*r + b*x + (2D)*x
proc stagD2*(sd:StaggeredD; r:Field; g:openArray[Field2];
             x:Field; a:SomeNumber; b:SomeNumber2) =
  template sf0:untyped = sd.sf
  template sb0:untyped = sd.sb
  let nd = g.len
  tic()
  for mu in 0..<nd:
    startSB(sf0[mu], x[ix])
  toc("startShiftF")
  for mu in 0..<nd:
    startSB(sb0[mu], g[mu][ix].adj*x[ix])
  toc("startShiftB")
  for ir in r[sd.subset]:
  #let ns0 = sd.subset.lowOuter
  #let ns1 = sd.subset.highOuter
  #let ns = ns1 - ns0
  #tFor iri, 0, <ns:
  #  let ir = ns0 + iri
    var rir{.noInit.}:type(r[ir])
    rir := a*r[ir] + b*x[ir]
    for mu in 0..<nd:
      localSB(sf0[mu], ir, imadd(rir, g[mu][ir], it), x[ix])
    for mu in 0..<nd:
      localSB(sb0[mu], ir, isub(rir, it), g[mu][ix].adj*x[ix])
      #var t{.noInit.}:type(load1(x[0]))
      #localSB(sb0[mu], ir, isub(rir, it), (mul(t,g[mu][ix].adj,x[ix]);t))
    assign(r[ir], rir)
  toc("local", flops=(18+nd*(72+66+6))*sd.subset.len)
  for mu in 0..<nd:
    boundarySB(sf0[mu], imadd(r[ir], g[mu][ir], it))
  toc("boundaryF")
  for mu in 0..<nd:
    boundarySB(sb0[mu], isub(r[ir], it))
  #threadBarrier()
  toc("boundaryB")

# r = m*x + sc*D*x
proc stagDN*(sd:openArray[StaggeredD]; r:openArray[Field]; g:openArray[Field2];
             x:openArray[Field]; m:SomeNumber; sc:SomeNumber=1.0) =
  stagDPN(sd, r, g, x, 6):
    #for i in 0..<n:
    rir := m*x[i][ir]
  #r[sd.subset] := (0.5*sc)*r

# r = m*x + sc*D*x
proc stagD*(sd:StaggeredD; r:Field; g:openArray[Field2];
            x:Field; m:SomeNumber; sc:SomeNumber=1.0) =
  stagD2(sd, r, g, x, 0, m/(0.5*sc))
  r[sd.subset] := (0.5*sc)*r
  #stagDP2(sd, r, g, x, 6):
  #  #for i in 0..<n:
  #  rir := m*getVec(x[ir], ic)

# r = m*x + sc*D*x
proc stagDb*(sd:StaggeredD; r:Field; g:openArray[Field2];
            x:Field; m:SomeNumber; sc:SomeNumber=1.0) =
  stagD2(sd, r, g, x, 0, m/(0.5*sc))
  #r[sd.subset] := (0.5*sc)*r
  #stagDP2(sd, r, g, x, 6):
  #  #for i in 0..<n:
  #  rir := m*getVec(x[ir], ic)

# r = m2 - Deo * Doe
proc stagD2ee*(sde,sdo:StaggeredD; r:Field; g:openArray[Field2];
               x:Field; m2:SomeNumber) =
  tic()
  var t{.global.}:type(x)
  if t==nil:
    threadBarrier()
    if threadNum==0:
      t = newOneOf(x)
    threadBarrier()
  #threadBarrier()
  #stagD(sdo, t, g, x, 0.0)
  toc("stagD2ee init")
  block:
    stagDP(sdo, t, g, x, 0):
      rir := 0
  toc("stagD2ee DP")
  threadBarrier()
  toc("stagD2ee barrier")
  #stagD(sde, r, g, t, 0.0)
  block:
    stagDM(sde, r, g, t, 6):
      rir := (4.0*m2)*x[ir]
  toc("stagD2ee DM")
  #threadBarrier()
  #r[sde.sub] := m2*x - r
  #for ir in r[sde.subset]:
  #  msubVSVV(r[ir], m2, x[ir], r[ir])
  #r[sde.sub] := 0.25*r

#[
# r = m2 - Deo * Doe
proc stagD2eeN*(sde,sdo:StaggeredD; r:Field; g:openArray[Field2];
                x:Field; m2:SomeNumber) =
  block:
    stagDPN(sdo, t, g, x, 0):
      rir := 0
  threadBarrier()
  block:
    stagDMN(sde, r, g, t, 6):
      rir := (4.0*m2)*x[ir]
]#

proc setBC*(g:openArray[Field]) =
  let gt = g[3]
  tfor i, 0..<gt.l.nSites:
    #let e = i div gt.l.nSitesInner
    if gt.l.coords[3][i] == gt.l.physGeom[3]-1:
      gt{i} *= -1
      #echoAll isMatrix(gt{i})
      #echoAll i, " ", gt[e][0,0]
proc stagPhase*(g:openArray[Field]) =
  const phases = [8,9,11,0]
  let l = g[0].l
  for mu in 0..<4:
    tfor i, 0..<l.nSites:
      var s = 0
      for k in 0..<4:
        s += (phases[mu] shr k) and l.coords[k][i].int
      if (s and 1)==1:
        g[mu]{i} *= -1
        #echoAll i, " ", gt[e][0,0]

proc newStag*[G,T](g:openArray[G];v:T):auto =
  var l = g[0].l
  template t:untyped =
    type(v[0])
  var r:Staggered[G,t]
  r.se = initStagDT(l, t, "even")
  r.so = initStagDT(l, t, "odd")
  r.g = @g
  r

proc newStag*[G](g:openArray[G]):auto =
  var l = g[0].l
  template t:untyped =
    type(l.ColorVector()[0])
    #SColorVectorV
  var r:Staggered[G,t]
  r.se = initStagDT(l, t, "even")
  r.so = initStagDT(l, t, "odd")
  r.g = @g
  r

proc newStag3*[G](g:openArray[G]):auto =
  var l = g[0].l
  template t:untyped =
    type(l.ColorVector()[0])
  var r:Staggered[G,t]
  r.se = initStagD3T(l, t, "even")
  r.so = initStagD3T(l, t, "odd")
  r.g = @g
  r

proc D*(s:Staggered; r,x:Field; m:SomeNumber) =
  stagD(s.se, r, s.g, x, m)
  stagD(s.so, r, s.g, x, m)
proc Ddag*(s:Staggered; r,x:Field; m:SomeNumber) =
  stagD(s.se, r, s.g, x, m, -1)
  stagD(s.so, r, s.g, x, m, -1)
proc eoReduce*(s:Staggered; r,b:Field; m:SomeNumber) =
  # r.even = (D^+ b).even
  #dump: "b.even.norm2"
  #dump: "b.odd.norm2"
  stagD(s.se, r, s.g, b, m, -1)
  #dump: r.even.norm2
  #dump: r.odd.norm2
proc eoReconstruct*(s:Staggered; r,b:Field; m:SomeNumber) =
  # r.odd = (b.odd - Doe r.even)/m
  stagD(s.so, r, s.g, r, 0.0, -1.0/m)
  r.odd += b/m

proc initSolverParams*():SolverParams =
  result.r2req = 1e-6
  result.maxits = 2000
  result.verbosity = 1
  result.subsetName = "even"

proc solve*(s:Staggered; r,x:Field; m:SomeNumber; sp0:SolverParams) =
  var sp = sp0
  sp.subset.layoutSubset(r.l, sp.subsetName)
  var t = newOneOf(r)
  var top = 0.0
  proc op(a,b:Field) =
    threadBarrier()
    if threadNum==0: top -= epochTime()
    stagD2ee(s.se, s.so, a, s.g, b, m*m)
    if threadNum==0: top += epochTime()
    #threadBarrier()
  threads:
    #echo "x2: ", x.norm2
    s.eoReduce(t, x, m)
    #echo "te2: ", t.even.norm2
  let t0 = epochTime()
  cgSolve(r, t, op, sp)
  let t1 = epochTime()
  threads:
    r[s.se.sub] := 4*r
    threadBarrier()
    s.eoReconstruct(r, x, m)
  let secs = t1-t0
  let flops = (s.g.len*4*72+60)*r.l.nEven*sp.finalIterations
  echo top
  echo "solve time: ", secs, "  Gflops: ", 1e-9*flops.float/secs
proc solve*(s:Staggered; r,x:Field; m:SomeNumber; res:float) =
  var sp = initSolverParams()
  sp.r2req = res
  #sp.maxits = 1000
  sp.verbosity = 1
  solve(s, r, x, m, sp)

proc solve2*(s:Staggered; r,x:Field; m:SomeNumber; res:float) =
  var sp:SolverParams
  sp.r2req = res
  sp.maxits = 100
  sp.verbosity = 1
  sp.subset.layoutSubset(r.l, "all")
  var t = newOneOf(r)
  proc op(a,b:Field) =
    #stagD2ee(s.se, s.so, r, s.g, x, m*m)
    threadBarrier()
    s.Ddag(t, b, m)
    #threadBarrier()
    #echo t.norm2
    s.D(a, t, m)
    #a := b
    #threadBarrier()
  #echo r.norm2
  cgSolve(r, x, op, sp)
  threads:
    s.Ddag(t, r, m)
    r := t

template foldl*(f,n,op:untyped):untyped =
  var r:type(f(0))
  r = f(0)
  for i in 1..<n:
    let
      a {.inject.} = r
      b {.inject.} = f(i)
    r = op
  r

when isMainModule:
  proc runtest(v1,v2,sdAll,sdEven,sdOdd,s,m:any) =
    let g = s.g
    let lo = g[0].l
    const nv = nVecs(v1[0])
    threads:
      v1 := 0
      #v2 := 1
      if myRank==0 and threadNum==0:
        when type(v1{0}) is AsVarVector:
          v1{0}[0] := 1
        else:
          v1{0} := 1
      threadBarrier()
      echo v1.norm2

      stagDb(sdAll, v2, g, v1, m)
      threadBarrier()
      echo v2.norm2
      #echo v2
      s.D(v2, v1, m)
      threadBarrier()
      echo v2.norm2

      for e in v1:
        template x(d:int):untyped = lo.vcoords(d,e)
        when type(v1{0}) is AsVarVector:
          v1[e][0].re := foldl(x, 4, a*10+b)
        else:
          for i in 0..<v1[e].ncols:
            v1[e][0,i].re := foldl(x, 4, a*10+b)
        #echo v1[e][0]
      threadBarrier()
      stagD(sdAll, v2, g, v1, 0.5)
      echo v1[0][0]
      echo v2[0][0]

    #let nrep = int(1e7/lo.physVol.float)
    #let nrep = int(2e8/lo.physVol.float)
    let nrep = int(1e9/lo.physVol.float)
    #let nrep = 1
    template makeBench(name:untyped; bar:untyped):untyped =
      proc `name T`(sd,v1,v2:any, ss="all") =
        resetTimers()
        var t0 = epochTime()
        threads:
          for rep in 1..nrep:
            stagDb(sd, v2, g, v1, 0.5)
            when bar: threadBarrier()
        var t1 = epochTime()
        let dt = t1-t0
        #var vol = lo.physVol.float
        var vol = lo.nSites.float
        if sd.sub != "all": vol *= 0.5
        let flops = nv * (6.0+g.len*2.0*72.0) * vol
        echo ss & "secs: ", dt, "  mf: ", (nrep.float*flops)/(1e6*dt)
        echoTimers()
      template name(sd:any, ss="all") = `name T`(sd, v1, v2, ss)
    subst(bench,_,benchB,_):
      makeBench(bench, false)
      makeBench(benchB, true)
      bench(sdAll, "all  ")
      benchB(sdAll, "all  ")
      bench(sdEven, "even ")
      benchB(sdEven, "even ")
      bench(sdOdd, "odd  ")
      benchB(sdOdd, "odd  ")
    proc benchEO() =
      resetTimers()
      var t0 = epochTime()
      threads:
        for rep in 1..nrep:
          stagD2ee(sdEven, sdOdd, v2, g, v1, 0.1)
      var t1 = epochTime()
      let dt = t1-t0
      #var vol = 0.5 * lo.physVol.float
      var vol = 0.5 * lo.nSites.float
      let flops = nv * (6.0+g.len*2.0*2.0*72.0) * vol
      echo "EO   secs: ", dt, "  mf: ", (nrep.float*flops)/(1e6*dt)
      #echoTimers()
    benchEO()

  qexInit()
  echo "rank ", myRank, "/", nRanks
  let cp = commandLineParams()
  #var lat = [4,4,4,4]
  var lat = [8,8,8,8]
  #var lat = [8,8,8,16]
  #var lat = [8,8,16,16]
  #var lat = [16,16,16,8]
  #var lat = [16,16,16,16]
  #var lat = [16,16,16,32]
  if cp.len>0:
    var i0 = 0
    if cp[0][0] notin {'0'..'9'}: inc i0
    for i in 0..<lat.len:
      lat[i] = (if (i0+i)<cp.len: parseInt(cp[i0+i]) else: lat[i-1])
  var lo = newLayout(lat)
  var v1 = lo.ColorVector()
  var v2 = lo.ColorVector()
  var g:array[4,type(lo.ColorMatrix())]
  for i in 0..<4:
    g[i] = lo.ColorMatrix()
    threads:
      g[i] := 1
  threads:
    g.setBC
    threadBarrier()
    for i in 0..<4:
      echo g[i].norm2
    threadBarrier()
    g.stagPhase
    threadBarrier()
    for i in 0..<4:
      echo g[i].norm2

  #g.loadGauge("l88.scidac")
  var sdAll = initStagD(v1, "all")
  var sdEven = initStagD(v1, "even")
  var sdOdd = initStagD(v1, "odd")
  var s = newStag(@g)
  var m = 0.1
  echo "done newStag"

  runtest(v1, v2, sdAll, sdEven, sdOdd, s, m)
  echoTimers()

  var sdAll3 = initStagD3(v1, "all")
  var sdEven3 = initStagD3(v1, "even")
  var sdOdd3 = initStagD3(v1, "odd")
  var g3:array[8,type(lo.ColorMatrix())]
  for i in 0..3:
    g3[2*i  ] = g[i]
    #g3[2*i+1] = g[i]
    g3[2*i+1] = lo.ColorMatrix()
    g3[2*i+1].random
  var s3 = newStag3(@g3)

  runtest(v1, v2, sdAll3, sdEven3, sdOdd3, s3, m)

  #[
  const nc = v1[0].len
  const nr = 8
  type MX* = Field[VLEN,MatrixArray[nr,nc,SComplexV]]
  #type MX* = Field[VLEN,MatrixArray[nr,nc,DComplexV]]
  var m1,m2:MX
  m1.new(lo)
  m2.new(lo)
  var sdAllM = initStagD(m1, "all")
  var sdEvenM = initStagD(m1, "even")
  var sdOddM = initStagD(m1, "odd")
  var sM = newStag(@g,m1)
  echo "testing multi matrix: ", nr
  stagD(sdAllM, m2, g, m1, m)
  runtest(m1, m2, sdAllM, sdEvenM, sdOddM, sM, m)
  echoTimers()
  ]#

  #[
  #const n = 4
  var n = 4
  if cp[0][0]=='n':
    n = parseInt(cp[0][1..^0])
  echo "n: ", n
  var v1a = newSeq[type(v1)](n)
  var v2a = newSeq[type(v2)](n)
  var sda = newSeq[type(sdAll)](n)
  var sda3 = newSeq[type(sdAll3)](n)
  #var sa = array[n,type(s)]
  v1a[0] = v1
  v2a[0] = v2
  sda[0] = sdAll
  sda3[0] = sdAll3
  #sda[0] = sdEven
  #sa[0] = s
  for i in 1..<n:
    v1a[i] = lo.ColorVector()
    v1a[i] := 1
    v2a[i] = lo.ColorVector()
    sda[i] = initStagD(v1, "all")
    sda3[i] = initStagD3(v1, "all")
    #sa[i] = newStag(@g)

  let nrep = int(2e8/lo.physVol.float)
  template makeBenchN(name:untyped; bar:bool):untyped =
    proc name(sd,g:any, ss="all") =
      resetTimers()
      var t0 = epochTime()
      threads:
        for rep in 1..nrep:
          stagDN(sd, v2a, g, v1a, 0.5)
          when bar: threadBarrier()
      var t1 = epochTime()
      let dt = t1-t0
      #var vol = lo.physVol.float
      var vol = lo.nSites.float
      if sd[0].sub != "all": vol *= 0.5
      let flops = n*(6.0+g.len*2.0*72.0) * vol
      echo ss & "secs: ", dt, "  mf: ", (nrep.float*flops)/(1e6*dt)
      #echoTimers()

  makeBenchN(benchN, false)

  stagDN(sda, v2a, g, v1a, m)
  benchN(sda, g)
  #echoTimers()

  #stagDN(sda3, v2a, g3, v1a, m)
  #benchN(sda3, g3)
  ]#

  qexFinalize()
