import ../mcmcTypes

import layout
import physics/[qcdTypes, stagSolve]

import strutils
import json

const
  ActionCGTol = 1e-20
  ForceCGTol = 1e-12
  ActionMaxCGIter = 10000
  ForceMaxCGIter = 10000

# Procs for creation of different kinds of staggered matter fields

proc checkJSON(info: JsonNode): JsonNode = 
  result = parseJson("{}")
  for key, keyVal in info: result[key] = keyVal
  if not result.hasKey("cg-tolerance-action"): result["cg-tolerance-action"] = %* ActionCGTol
  if not result.hasKey("cg-maxits-action"): result["cg-maxits-action"] = %* ActionMaxCGIter
  if not result.hasKey("cg-tolerance-force"): result["cg-tolerance-force"] = %* ForceCGTol
  if not result.hasKey("cg-maxits-force"): result["cg-maxits-force"] = %* ForceMaxCGIter
  result["field-type"] = %* "staggered"

proc newStaggeredField(self: var LatticeField; info: JsonNode) = 
  self.stagActionSolverParams = initSolverParams()
  self.stagForceSolverParams = initSolverParams()
  self.stagActionSolverParams.r2req = info["cg-tolerance-action"].getFloat()
  self.stagActionSolverParams.maxits = info["cg-maxits-action"].getInt()
  self.stagForceSolverParams.r2req = info["cg-tolerance-force"].getFloat()
  self.stagForceSolverParams.maxits = info["cg-maxits-force"].getInt()
  self.staggeredMasses = newSeq[float]()

proc newStaggeredFermion*(l: Layout; staggeredInformation: JsonNode): auto = 
  var 
    info = checkJSON(staggeredInformation)
    stream = newMCStream("new staggered fermion")
  if not info.hasKey("mass"): qexError "Must specify fermion mass"

  result = l.newLatticeField(info)
  result.newStaggeredField(info)
  result.staggeredFields = newSeq[l.TT]()

  result.staggeredAction = StaggeredFermion
  result.staggeredFields.add l.ColorVector()
  result.staggeredMasses.add info["mass"].getFloat()

  stream.add "  mass = " & $(info["mass"].getFloat())
  stream.finishStream

proc newStaggeredBoson*(l: Layout; staggeredInformation: JsonNode): auto = 
  var 
    info = checkJSON(staggeredInformation)
    stream = newMCStream("new staggered boson")
  if not info.hasKey("mass"): qexError "Must specify boson mass"

  result = l.newLatticeField(info)
  result.newStaggeredField(info)
  result.staggeredFields = newSeq[l.TT]()

  result.staggeredAction = StaggeredBoson
  result.staggeredFields.add l.ColorVector()
  result.staggeredMasses.add info["mass"].getFloat()

  stream.add "  mass = " & $(info["mass"].getFloat())
  stream.finishStream

proc newStaggeredHasenbuschFermion*(
    l: Layout; 
    staggeredInformation: JsonNode
  ): auto = 
  var 
    info = checkJSON(staggeredInformation)
    stream = newMCStream("new staggered Hasenbusch fermion")
  if not info.hasKey("mass1"): qexError "Must specify first Hasenbusch mass"
  if not info.hasKey("mass2"): qexError "Must specify second Hasenbusch mass"
  let masses = @[info["mass1"].getFloat(), info["mass2"].getFloat()]

  result = l.newLatticeField(info)
  result.newStaggeredField(info)
  result.staggeredFields = newSeq[l.TT]()

  result.staggeredAction = StaggeredFermion
  for mass in masses:
    result.staggeredFields.add l.ColorVector()
    result.staggeredMasses.add mass

  stream.add "  mass1 = " & $(info["mass1"].getFloat())
  stream.add "  mass2 = " & $(info["mass2"].getFloat())
  stream.finishStream

# Generic procs for operations w/ staggered matter fields

proc phi(self: LatticeField): auto = self.staggeredFields[0]

proc phi1(self: LatticeField): auto = self.staggeredFields[0]

proc phi2(self: LatticeField): auto = self.staggeredFields[1]

proc mass(self: LatticeField): float = self.staggeredMasses[0]

proc mass1(self: LatticeField): float = self.staggeredMasses[0]

proc mass2(self: LatticeField): float = self.staggeredMasses[1]

proc zero(phi: auto) =
  threads: phi := 0

proc zeroOdd(phi: auto) =
  threads: phi.odd := 0

proc sq(x: float): float = x*x

proc normSquared(psi: auto): float =
  var nrm2 = 0.0
  threads:
    let psi2 = psi.norm2
    threadBarrier()
    threadMaster: nrm2 = psi2
  result = nrm2

proc applyD(
    D: DiracOperator;
    psi: auto;
    phi: auto;
    mass: float;
  ) = 
  threads: D(D.stag, psi, phi, mass)

proc applyDdag(
    D: DiracOperator;
    psi: auto;
    phi: auto;
    mass: float;
  ) =
  threads: D(D.stag, psi, phi, -mass)

proc applyNegDdagOdd(
    D: DiracOperator;
    psi: auto;
    phi: auto;
  ) =
  threads:
    stagD2(D.stag.so, psi, D.stag.g, phi, 0, 0)
    threadBarrier()
    psi.odd := -0.5*psi
    psi.even := 0

proc applyDdag2OddAndReplaceEven(
    D: DiracOperator;
    psi: auto;
    phi: auto;
  ) =
  threads:
    stagD2(D.stag.so, psi, D.stag.g, phi, 0, 0)
    threadBarrier()
    psi.even := phi

proc solveD*(
    D: DiracOperator;
    psi: auto;
    phi: auto;
    mass: float;
    sp0: var SolverParams
  ) =
  psi.zero
  if mass > Small64: solve(D.stag, psi, phi, mass, sp0)
  else:
    var sp = sp0
    sp.resetStats()
    sp.verbosity = sp0.verbosity
    sp.usePrevSoln = false

    threads: psi := 0
    solveEE(D.stag, psi, phi, 0, sp)
    threads: psi.even := 4*psi

    sp.calls = 1
    sp0.addStats(sp)

proc solveDdag*(
    D: DiracOperator;
    psi: auto;
    phi: auto;
    mass: float;
    sp0: var SolverParams
  ) =
  psi.zero
  if mass > epsilon(float): solve(D.stag, psi, phi, -mass, sp0)
  else:
    var sp = sp0
    sp.resetStats()
    sp.verbosity = sp0.verbosity
    sp.usePrevSoln = false

    threads: psi := 0
    solveEE(D.stag, psi, phi, 0, sp)
    threads: psi.even := 4*psi

    sp.calls = 1
    sp0.addStats(sp)

proc outer(f: auto; psi: auto; shifter: auto; dtau: float) =
    let n = psi[0].len
    threads:
      for mu in 0..<f.len:
        for s in f[mu]:
          forO a, 0, n-1:
            forO b, 0, n-1:
              f[mu][s][a,b] += dtau * psi[s][a] * shifter[mu].field[s][b].adj

#[ 
Methods for 
1.) getting phi, 
2.) calculating action, 
3.) calculating (partial) force
]#

proc getStaggeredField*(self: var LatticeField; D: DiracOperator) =
  case self.staggeredAction:
    of StaggeredFermion,StaggeredBoson: zero(self.phi)
    of StaggeredHasenbuschFermion:
      zero(self.phi1)
      zero(self.phi2)
  case self.staggeredAction:
    of StaggeredFermion: D.applyDdag(self.phi, D.stagPsi, self.mass)
    of StaggeredHasenbuschFermion: 
      D.solveDdag(self.phi2, D.stagPsi, self.mass2, self.stagActionSolverParams)
      D.applyDdag(self.phi1, self.phi2, self.mass1)
    of StaggeredBoson: 
      D.solveDdag(self.phi, D.stagPsi, self.mass, self.stagActionSolverParams)
  zeroOdd(self.phi)

proc staggeredAction*(self: var LatticeField; D: var DiracOperator): float =
  D.stagPsi.zero
  case self.staggeredAction:
    of StaggeredFermion: 
      D.solveDdag(D.stagPsi, self.phi, self.mass, self.stagActionSolverParams)
      if self.mass <= Small64: D.applyNegDdagOdd(D.stagPsi, D.stagPsi)
    of StaggeredHasenbuschFermion:
      zero(self.phi2)
      D.applyDdag(self.phi2, self.phi1, self.mass2)
      D.solveDdag(D.stagPsi, self.phi2, self.mass1, self.stagActionSolverParams)
    of StaggeredBoson: 
      D.applyDdag(D.stagPsi, self.phi, self.mass)
  result = 0.5*D.stagPsi.normSquared

proc kineticAction(p: auto): float =
  var p2: float
  threads:
    var p2r = 0.0
    for mu in 0..<p.len: p2r += p[mu].norm2
    threadBarrier()
    threadMaster: p2 = p2r
  result = 0.5*p2 - 16.0*p[0].l.physVol

proc stagForce*(
    self: var LatticeField;
    D: var DiracOperator;
    fdtau: float;
    f: auto
  ) = 
  var dtau = -0.5 * fdtau
  D.stagPsi.zero
  case self.staggeredAction:
    of StaggeredFermion:
      D.solveD(D.stagPsi, self.phi, self.mass, self.stagForceSolverParams)
      if self.mass < Small64: 
        D.applyDdag2OddAndReplaceEven(D.stagPsi, D.stagPsi)
    of StaggeredHasenbuschFermion: 
      D.solveD(D.stagPsi, self.phi1, self.mass1, self.stagForceSolverParams)
    of StaggeredBoson: 
      D.applyDdag2OddAndReplaceEven(D.stagPsi, self.phi)
  case self.staggeredAction:
    of StaggeredFermion, StaggeredHasenbuschFermion, StaggeredBoson:
      for mu in 0..<f.len: discard D.stagShifter[mu] ^* D.stagPsi
  case self.staggeredAction:
    of StaggeredFermion:
      if self.mass > Small64: dtau = dtau / self.mass
      else: dtau = -0.5 * dtau
    of StaggeredHasenbuschFermion: 
      dtau = dtau * (self.mass2.sq - self.mass1.sq) / self.mass1
    of StaggeredBoson: dtau = 0.5 * dtau
  case self.staggeredAction:
    of StaggeredFermion, StaggeredHasenbuschFermion, StaggeredBoson:
      f.outer(D.stagPsi, D.stagShifter, dtau)
  