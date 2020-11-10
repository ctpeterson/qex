import qex
import gauge, physics/qcdTypes
import os, strutils, times

proc sumEnergy(fr,fi:any, J, h:any, g:any, p:seq[float], sf,sb:any) =
  fr := 0
  fi := 0
  for nu in 0..<g.l.nDim:
    discard sf[nu] ^* g
    discard sb[nu] ^* g
    var pf,cpf,spf,pb,cpb,spb:typeof(g[0])
    threadBarrier()
    for i in fr:
      pf := sf[nu].field[i] - p[nu]
      cpf = cos pf
      spf = sin pf
      pb := sb[nu].field[i] + p[nu]
      cpb = cos pb
      spb = sin pb
      fr[i] += cpf + cpb
      fi[i] += spf + spb
  fr := J*fr + h
  fi := J*fi

type PhaseDiff[F,E] = object
  cosd,sind:seq[float]
  f: seq[Shifter[F,E]]

proc phaseDiff(del:var PhaseDiff,g:any,p:seq[float],tdir:seq[bool]):auto =
  let
    # del cannot be captured by nim in threads
    f = del.f
    cosd = cast[ptr UnCheckedArray[float]](del.cosd[0].addr)
    sind = cast[ptr UnCheckedArray[float]](del.sind[0].addr)
  threads:
    for nu in 0..<g.l.nDim:
      if not tdir[nu]: continue
      var d,t,s:typeof(g[0])
      discard f[nu] ^* g
      threadBarrier()
      for i in g:
        d := f[nu].field[i] - g[i] - p[nu]
        t += cos(d)
        s += sin(d)
      var
        v = t.simdSum
        u = s.simdSum
      v.threadRankSum
      u.threadRankSum
      threadSingle:
        cosd[nu] = v
        sind[nu] = u

type HeatBath[F,E] = object
  fr,fi: F
  sf,sb: array[2,seq[Shifter[F,E]]]
  subs: array[2,Subset]
  del: PhaseDiff[F,E]

proc newHeatBath(lo:any):auto =
  let
    nd = lo.nDim
    fr = lo.Real
    fi = lo.Real
  type
    F = typeof(fr)
    E = typeof(fr[0])
  var r = HeatBath[F,E](fr:fr, fi:fi)
  const p = ["even","odd"]
  for j in 0..1:
    r.sf[j] = newseq[Shifter[F,E]](nd)
    r.sb[j] = newseq[Shifter[F,E]](nd)
    for i in 0..<nd:
      r.sf[j][i] = newShifter(fr, i, 1, p[j])
      r.sb[j][i] = newShifter(fr, i, -1, p[j])
  r.subs[0].layoutSubset(lo,"e")
  r.subs[1].layoutSubset(lo,"o")
  r.del.cosd = newseq[float](nd)
  r.del.sind = newseq[float](nd)
  r.del.f = newseq[Shifter[F,E]](nd)
  for i in 0..<nd:
    r.del.f[i] = newShifter(fr, i, 1)
  r

proc evolve(H:HeatBath, g:any, d:var seq[float], tdir:seq[bool], gc:any, r:any, R:var RngMilc6,
    sample = true, twistSample = true, jump = true, twistJump = true) =
  tic("heatbath")
  let
    lo = g.l
    nd = lo.nDim
    (beta, J, h) = gc
    p = d
    z = newseq[float](nd)
  if H.subs.len != 2:
    qexError "HeatBath only works with even-odd subsets for now."
  if sample:
    tic("threads")
    threads:
      for j in 0..<H.subs.len:
        let
          s = H.subs[j]
          so = H.subs[(j+1) mod 2]
        sumEnergy(H.fr[s], H.fi[s], J, h, g, p, H.sf[j], H.sb[j])
        threadBarrier()
        for i in g[s].sites:
          let
            yr = H.fr{i}[][]
            yi = H.fi{i}[][]
            lambda = beta*hypot(yi, yr)
            phi = arctan2(yi, yr)
          g{i} := vonMises(r{i}, lambda)+phi
    toc("sample")
  if twistSample:
    tic()
    var del = H.del
    del.phaseDiff(g,z,tdir)
    for mu in 0..<nd:
      if not tdir[mu]: continue
      let
        yr = J*del.cosd[mu] + h
        yi = J*del.sind[mu]
        phi = arctan2(yi,yr)
      d[mu] = vonMises(R, beta*hypot(yi,yr))+phi
    toc("twist sample")
  if jump:
    tic("threads")
    threads:
      # over relaxation: flip
      for j in 0..<H.subs.len:
        let
          s = H.subs[j]
          so = H.subs[(j+1) mod 2]
        sumEnergy(H.fr[s], H.fi[s], J, h, g, p, H.sf[j], H.sb[j])
        threadBarrier()
        for i in g[s].sites:
          let
            yr = H.fr{i}[][]
            yi = H.fi{i}[][]
          g{i} := 2.0*arctan2(yi,yr)-g{i}[][]
    toc("flip")
  if twistJump:
    tic()
    var del = H.del
    del.phaseDiff(g,z,tdir)
    for mu in 0..<nd:
      if not tdir[mu]: continue
      let
        yr = J*del.cosd[mu] + h
        yi = J*del.sind[mu]
        phi = arctan2(yi,yr)
      d[mu] = 2.0*phi-d[mu]
    toc("twist flip")
  toc("end")

proc magnet(g:any):auto =
  tic("magnet")
  var mr,mi = 0.0
  threads:
    var t,s:typeof(g[0])
    for i in g:
      t += cos(g[i])
      s += sin(g[i])
    var
      v = t.simdSum
      u = s.simdSum
    v.threadRankSum
    u.threadRankSum
    threadSingle:
      mr = v
      mi = u
  toc("done")
  (mr, mi)

proc showMeasure[F,E](del:var PhaseDiff[F,E],g:F,d:seq[float],label="") =
  let
    (mr,mi) = g.magnet
    v = 1.0/g.l.physVol.float
    s = (mr*mr+mi*mi)*v
    nd = g.l.nDim
  var tdir = newseq[bool](nd)
  for i in 0..<nd: tdir[i] = true
  del.phaseDiff(g,d,tdir)
  echo label,"magnet: ",mr," ",mi," ",s
  var diff = ""
  for i in 0..<nd:
    diff &= "CosSinDel" & $i & ": " & $(del.cosd[i]*v) & " " & $(del.sind[i]*v) & "\t"
  diff.setlen(diff.len-1)
  echo label,diff

proc showTwist[T](d:openarray[T],label="") =
  var pd = label & "twist: "
  for mu in 0..<d.len: pd.add " " & $d[mu]
  echo pd

proc hitFreq(num, freq:int):bool = freq>0 and 0==num mod freq

proc toSeqBool(s:seq[int]):seq[bool] =
  result = newseq[bool](s.len)
  for i in 0..<s.len:
    result[i] = if s[i]==0: false else: true

qexinit()
tic()

letParam:
  #lat = @[8,8,8,8]
  #lat = @[8,8,8]
  lat = @[32,32]
  #lat = @[1024,1024]
  beta = 1.0
  J = 1.0
  h = 0.0
  sweeps = 10
  sampleFreq = 1
  jumpFreq = 1
  twistSampleFreq = 1
  twistJumpFreq = 1
  twistDirs:toSeqBool = @[1,0]
  twistAngle = @[0.0,0.0]
  measureFreq = 1
  seed:uint64 = int(1000*epochTime())

echoParams()
echo "rank ", myRank, "/", nRanks
threads: echo "thread ", threadNum, "/", numThreads

if twistAngle.len != lat.len:
  qexError "The lengths of twistAngle(" & $twistAngle & ") and lat(" & $lat & ") differ."
if twistAngle.len != twistDirs.len:
  qexError "The lengths of twistAngle(" & $twistAngle & ") and twistDirs(" & $twistDirs & ") differ."

let
  gc = (beta:beta, J:J, h:h)
  lo = lat.newLayout
  vol = lo.physVol

var
  r = lo.newRNGField(RngMilc6, seed)
  R:RngMilc6  # global RNG for the twisting angle
  g = lo.Real
  d = twistAngle
  H = lo.newHeatBath
R.seed(seed,987654321)

threads:
  for i in g.sites:
    let u = uniform r{i}
    g{i} := PI*(2.0*u-1.0)

d.showTwist("Initial: ")
H.del.showMeasure(g,d, "Initial: ")

toc("init")

for n in 1..sweeps:
  tic("sweep")
  echo "Begin sweep: ",n

  H.evolve(g,d,twistDirs,gc,r,R,
    hitFreq(n,sampleFreq),
    hitFreq(n,twistSampleFreq),
    hitFreq(n,jumpFreq),
    hitFreq(n,twistJumpFreq))
  toc("evolve")
  for mu in 0..<lo.nDim:
    if twistDirs[mu]:
      d[mu] = floormod(d[mu]+PI,2*PI) - PI

  d.showTwist
  if hitFreq(n,measureFreq): H.del.showMeasure(g,d)
  toc("measure")

toc("done")
echoTimers()
qexfinalize()
