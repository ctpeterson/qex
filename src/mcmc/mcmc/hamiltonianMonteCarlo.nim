import ../mcmcTypes
import ../utilities/gaugeUtils
import ../actions/latticeAction
import ../actions/latticeActionUtils
import randomNumberGeneration
import molecularDynamics

import options
import math

export latticeAction

proc kineticAction(p: auto): float =
  var p2: float
  threads:
    var p2r = 0.0
    for mu in 0..<p.len: p2r += p[mu].norm2
    threadBarrier()
    threadMaster: p2 = p2r
  result = 0.5*p2 - 16.0*p[0].l.physVol

template runHamiltonianMonteCarlo*(
    self: var LatticeFieldTheory; 
    samples: int;
    measurements: untyped 
  ) =

  proc hamiltonian: float {.gensym.} =
    result = self.p.kineticAction
    for action in self.actions: 
      result += action.getAction

  proc trajectory {.gensym.} =
    self.actions.trajectory(0.0, self.u[], self.f, self.p, type(self.actions[0]))

  proc metropolis(hi,hf: float): bool {.gensym.} =
    result = false
    case self.sRNG.uniform <= exp(hi-hf):
      of true:
        reunit(self.u[])
        result = true
      of false: self.u[].setGauge(self.bu)

  template u: untyped {.inject.} = self.u[]

  case self.start:
    of UnitGauge: unit(u)
    of RandomGauge: self.pRNG.warm(u)
    of ReadGauge: discard

  for smp in 0..<samples:
    var
      hi {.inject.} = 0.0
      hf {.inject.} = 0.0

    # Back up gauge field
    self.bu.setGauge(self.u[])

    # Smear gauge fields
    for action in self.actions: action.smearGauge

    # Momentum heatbath
    self.pRNG.randomTAHGaussian(self.p)

    # Fermion heatbath
    for action in self.actions: 
      action.fermionHeatbath

    # Get Hi
    hi = hamiltonian()

    # Evolve
    trajectory()

    # Smear gauge fields
    for action in self.actions: action.smearGauge

    # Get Hf
    hf = hamiltonian()

    # Metropolis
    let accepted {.inject.} = metropolis(hi,hf)

    # Let user make their own measurements
    let sample {.inject,used.} = smp
    measurements
  
if isMainModule:
  qexInit()

  var
    ftInfo = %* {
      "lattice-geometry": @[4,4,4,4],
      "mpi-geometry": @[1,1,1,1],
      "monte-carlo-algorithm": "hmc",
      "trajectory-length": 1.0,
      "serial-random-number-seed": 987654321,
      "parallel-random-number-seed": 987654321,
      "serial-random-number-generator": "milc",
      "parallel-random-number-generator": "milc",
      "start": "cold"
    }
    aInfo = %* {
      "smearing": "nhyp",
      "smearing-coefficients": @[0.4,0.5,0.5],
      "boundary-conditions": "aaaa"
    }
    gInfo = %* {
      "action": "Wilson",
      "beta": 9.0,
      "steps": 10,
      "integrator": "2MN"
    }
    fInfo = %* {
      "mass": 0.0,
      "integrator": "2MN",
      "steps": 10
    }
    bInfo = %* {
      "mass": 0.75,
      "integrator": "2MN",
      "steps": 10
    }
    sbInfo = %* {
      "mass": 0.75,
      "integrator": "2MN",
      "steps": 2
    }

  var hmc = newLatticeFieldTheory(ftInfo):
    fieldTheory.addGaugeAction(gInfo)
    fieldTheory.addMatterAction(aInfo):
      action.addStaggeredFermion(fInfo)
      action.addStaggeredBoson(bInfo)
  #      subAction.addStaggeredBoson(sbInfo) # Nested
  #      subAction.addStaggeredBoson(sbInfo) # Nested

  hmc.runHamiltonianMonteCarlo(2):
    echo sample, " ", accepted, " ", hf-hi
    plaquette(u)
    polyakov(u)

  qexFinalize()