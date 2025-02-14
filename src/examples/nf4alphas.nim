import qex
import examples/[hisqhmc_h]
import sequtils,parseutils,strutils
import parseopt,json

const banner = """
|---------------------------------------------------------------|
 Quantum EXpressions (QEX)

 QEX authors: James Osborn & Xiao-Yong Jin
 HISQ HMC authors: 
   - Curtis Taylor Peterson [C.T.P.] (Michigan State University)
   - James Osborn (Argonne National Laboratory)
   - Xiao-Yong Jin (Argonne National Laboratory) 
 QEX GitHub: https://github.com/jcosborn/qex
 C.T.P. email: curtistaylorpetersonwork@gmail.com
 cite: Proceedings of Science (PoS) LATTICE2016 (2017) 271
|---------------------------------------------------------------|
"""

qexInit()

echo banner
echo "rank ", myRank, "/", nRanks
threads: echo "thread ", threadNum, "/", numThreads

let
  prompt = readCMD()
  baseFilename = prompt["ensemble"].getStr()

# Proc for calculating plaquette
proc plaquette[T](u: T) =
  let
    pl = u.plaq
    nl = pl.len div 2
    ps = pl[0..<nl].sum * 2.0
    pt = pl[nl..^1].sum * 2.0
    ptot = 0.5*(ps+pt)
  echo "MEASplaq ss: ",ps,"  st: ",pt,"  tot: ",ptot
  
# Proc for calculating Polyakov loop
proc polyakov[T](u: T) =
  let pg = u[0].l.physGeom
  var pl = newseq[typeof(u.wline @[1])](pg.len)
  for i in 0..<pg.len: pl[i] = u.wline repeat(i+1, pg[i])
  let
    pls = pl[0..^2].sum / float(pl.len-1)
    plt = pl[^1]
  echo "MEASploop spatial: ",pls.re," ",pls.im," temporal: ",plt.re," ",plt.im
  
# Proc for measuring chiral condensate
proc condensate(hmc: auto) =
  var
    pbpsp: SolverParams
    tmpa = hmc.stag.g[0].l.ColorVector()
    tmpb = hmc.stag.g[0].l.ColorVector()
    pbptote = 0.0
    pbptoto = 0.0
  let 
    mass = hmc.mass
    vol = hmc.stag.g[0].l.physVol.float
    nsources = hmc.jsonInfo["measurements"]["chiral-condensate"]["sources"].getInt()
  pbpsp.r2req = ActionCGTol
  pbpsp.maxits = ActionMaxCGIter
  for source in 0..<nsources:
    threads: tmpa.gaussian(hmc.prng.milc)
    hmc.stag.solve(tmpb,tmpa,mass,pbpsp)
    threads:
      let 
        pbpe = 0.5*mass*tmpb.even.norm2/vol
        pbpo = 0.5*mass*tmpb.odd.norm2/vol
      threadMaster: echo "MEASpbp (",source,") mass: ",mass," pbpe: ",pbpe," pbpo: ",pbpo
      pbptote += pbpe/float(nsources)
      pbptoto += pbpo/float(nsources)
  echo "MEASpbp (avg) mass: ",mass," pbpe: ",pbptote," pbpo: ",pbptoto

# Construct HMC object
var hmc = newHisqHMC:
  # Gauge link update
  proc mdt(dtau: float) = hisq.updateGauge(dtau)

  # Momentum update
  proc mdvAll(dtau: openarray[float]) =
    let (dtauG,dtauF) = (dtau[0],dtau[1])
    if (dtauG != 0.0): hisq.updateMomentumGauge(dtauG)
    if (dtauF != 0.0): hisq.updateMomentumFermion(dtauF)

  # Construct integrator according to mdEvolve scheme
  let 
    (VAll,T) = newIntegratorPair(mdvAll,mdt)
    (V,Vf) = (VAll[0],VAll[1])
  integrator = newParallelEvolution(
    gaugeIntegrator(steps = gaugeSteps, V = V, T = T),
    fermionIntegrator(steps = fermionSteps, V = Vf, T = T)
  )

  # Read information from disk
  if start == "read":
    let fn = baseFilename & "_" & $(hisq.traj0)
    hisq.readGauge(fn & ".lat")
    if hisq.traj0 > 0:
      hisq.readSerialRNG(fn & ".serialRNG")
      hisq.readParallelRNG(fn & ".parallelRNG")
    u.plaquette
    u.polyakov
    u.reunit

# Do HMC
echo $(hmc)
hmc.sample:
  hmc.prepare()
  hmc.evolve()
  hmc.finish:
    let output = $(info.dH) & " exp(dH): " & $(info.expdH) & " rand: " & $(info.rnd)
    case accepted:
      of true: echo "ACC: ", output
      of false: echo "REJ: ", output
    if hmc.jsonInfo.hasKey("measurements"):
      if hmc.jsonInfo["measurements"].hasKey("plaquette"): u.plaquette
      if hmc.jsonInfo["measurements"].hasKey("polyakov"): u.polyakov
      if hmc.jsonInfo["measurements"].hasKey("chiral-condensate"): hmc.condensate
    if hmc.jsonInfo.hasKey("checkpoint"):
      let saveFreq = hmc.jsonInfo["checkpoint"]["frequency"].getInt()
      if (saveFreq > 0) and (((trajectory + 1) mod saveFreq) == 0):
        let fn = baseFilename & "_" & $(trajectory + 1)
        hmc.writeGauge(fn & ".lat")
        hmc.writeSerialRNG(fn & ".serialRNG")
        hmc.writeParallelRNG(fn & ".parallelRNG")

qexFinalize()