import qex
import layout
import gauge/[hypsmear, stoutsmear]

import json
import options

import ../mcmcTypes
import ../fields/gaugeFields
import ../fields/staggeredFields

# Procs/templates for constructing full action

proc checkJSON(l: Layout; info: JsonNode; actionType: ActionType): JsonNode = 
  result = parseJson("{}")
  for key, keyVal in info: result[key] = keyVal
  case actionType:
    of PureGauge: result["action-type"] = %* "gauge"
    of PureMatter: result["action-type"] = %* "matter"
    of GaugeMatter: result["action-type"] = %* "gauge-matter"
  case actionType:
    of PureGauge: discard
    of PureMatter, GaugeMatter:
      if not result.hasKey("smearing"): 
        result["smearing"] = %* "nhyp"
        result["smearing-coefficients"] = %* [0.4,0.5,0.5]
      if not result.hasKey("boundary-conditions"):
        var bcs = ""
        for dim in 0..<l.nDim-1: bcs = bcs & "p"
        bcs = bcs & "a"

template add(self: var LatticeAction; matter: untyped) =
  self.fields.add matter
  let field = self.fields[^1].field
  case field:
    of StaggeredMatterField, WilsonMatterField:
      if not self.D.hasKey(field): 
        self.D[field] = newDiracOperator(self.smear.su, field)
    of Gauge: discard

proc newLatticeAction*(
    l: Layout;
    actionInfo: JsonNode
  ): auto = 
  var info = l.checkJSON(actionInfo, PureGauge)
  result = l.newLatticeAction(info, l.SS, l.TT, l.UU, l.VV, l.WW, l.XX)
  result.add l.newGaugeAction(info)

template newLatticeAction*(
    l: Layout;
    actionInfo: JsonNode;
    matterFields: varargs[untyped]
  ): auto = 
  var info = l.checkJSON(actionInfo, PureMatter)
  result = l.newLatticeAction(info, l.SS, l.TT, l.UU, l.VV, l.WW, l.XX)
  for matter in matterFields: result.add matter

template newLatticeAction*(
    l: Layout;
    actionInfo, gaugeInfo: JsonNode; 
    matterFields: varargs[untyped]
  ): auto = 
  var info = l.checkJSON(actionInfo, GaugeMatter)
  result = l.newLatticeAction(info, l.SS, l.TT, l.UU, l.VV, l.WW, l.XX)
  result.add l.newGaugeAction(gaugeInfo)
  for matter in matterFields: result.add matter

# Helper procs for action/force calculation

proc rephase(u: auto) =
  threads: u.stagPhase

proc setBoundaryConditions(self: LatticeAction; u: auto) =
  let bcs = self.boundaryConditions
  threads:
    for mu in 0..<u.len:
      if $bcs[mu] == $"a":
        tfor i, 0..<u[mu].l.nSites:
          if u[mu].l.coords[mu][i] == u[mu].l.physGeom[mu]-1:
            u[mu]{i} *= -1.0

# Procs for calculating action & forces

proc action*(self: var LatticeAction): float = 
  result = self.fields["gauge"].gaugeAction()

proc action*[S](
    self: var LatticeAction; 
    u: seq[S]; 
    parallelRNG: Option[ParallelRNG] = none(ParallelRNG)
  ): Table[string,float] =
  result = initTable[string,float]()
  case self.action:
    of PureGauge: result["gauge"] = self.fields["gauge"].gaugeAction(u)
    of PureMatter, GaugeMatter:
      #[ !!!!!! CALCULATE CONTRIBUTION FROM GAUGE ACTION !!!!!! ]#

      # Switches for matter fields
      var rephased = false

      # Take care of smearing
      case self.smearing:
        of Hypercubic: 
          self.smear.nhyp.smear(u, self.smear.su, self.smear.nHYPInfo)
          self.setBoundaryConditions(self.smear.su)
        of Stout: discard
        of NoSmearing: self.setBoundaryConditions(u)

      # Calculate contribution from staggered matter field
      for key, _ in self.fields:
        if (self.fields[key].field == StaggeredMatterField):
          # Rephase
          if not rephased:
            case self.smearing:
              of Hypercubic, Stout: rephase(self.smear.su)
              of NoSmearing: rephase(u)
            rephased = true

          # Optional fermion heatbath
          if not parallelRNG.isNone:
            randomComplexGaussian(parallelRNG.get, self.stagPsi)
            self.fields[key].getStaggeredField(self.stagD, self.stagPsi)

          # Calculate & append contribution for full action
          result[key] = self.fields[key].staggeredAction(self.stagD, self.stagPsi)
      if (self.smearing == NoSmearing) and (rephased): rephase(u)

      #[Wilson matter field - not implemented]#

      # Undo boundary condition rephase if no smearing
      case self.smearing:
        of Hypercubic, Stout: discard
        of NoSmearing: self.setBoundaryConditions(u)

proc updateMomentum*[S](self: var LatticeAction; f,p: seq[S]) =
  # Switches for momentum update
  var updateMomentum = false

  # Momentum update from gauge field
  let dtau =  self.fields["gauge"].dtau[0]
  if dtau > epsilon(float): updateMomentum = true
  if updateMomentum: 
    self.fields["gauge"].gaugeForce(f)
    threads:
      for mu in 0..<f.len:
        for s in f[mu]: p[mu][s] -= dtau*f[mu][s]

proc updateMomentum*[S](self: var LatticeAction; u: seq[S]; f,p: seq[S]) =
  # Switches for momentum update
  var updateMomentum = false

  # Calculate force
  case self.action:
    of PureGauge:
      let dtau =  self.fields["gauge"].dtau[0]
      if dtau > epsilon(float): updateMomentum = true
      if updateMomentum: 
        self.fields["gauge"].gaugeForce(u,f)
        threads:
          for mu in 0..<f.len:
            for s in f[mu]: f[mu][s] := dtau*f[mu][s]
    of PureMatter:
      # Switches for matter fields
      var fieldsRephased = false

      # Check if smearing & momentum update need to be done
      for key, _ in self.fields:
        if self.fields[key].dtau[0] > epsilon(float): 
          updateMomentum = true

      # Calculate force
      if updateMomentum:
        let sft = StaggeredMatterField

        # Zero force
        threads:
          for mu in 0..<f.len: 
            for s in f[mu]: f[mu][s] := 0

        # Smear fields
        case self.smearing:
          of Hypercubic:
            self.smear.smearedForce = smearGetForce(
              self.smear.nhyp, u, self.smear.su, self.smear.nHYPInfo
            )
            self.setBoundaryConditions(self.smear.su)
          of Stout: discard
          of NoSmearing: self.setBoundaryConditions(u)

        # Take care of contribution to force from staggered fields
        for key, _ in self.fields:
          let dtau = self.fields[key].dtau[0]
          let field = self.fields[key].field
          if (dtau != 0.0) and (field == sft):
            if not fieldsRephased: # Rephase gauge fields - done only once 
              case self.smearing:
                of Hypercubic: rephase(self.smear.su)
                of Stout: discard
                of NoSmearing: rephase(u)
              fieldsRephased = true
            self.fields[key].staggeredPartialForce(self.D[sft], self.stagPsi, f)
        threads: # Rephase force
          f.stagPhase
          threadBarrier()
          for mu in 0..<f.len:
            for s in f[mu].odd: f[mu][s] *= -1
        case self.smearing: # Undo rephase if no smearing
          of Hypercubic, Stout: discard
          of NoSmearing:
            if fieldsRephased: rephase(u)

        #[Update from Wilson fields - not implemeted, will need to undo su rephase]#

        # Set boundary conditions, smear, and project
        self.setBoundaryConditions(f)
        case self.smearing: 
          of Hypercubic:
            self.smear.smearedForce(f,f)
            threads:
              for mu in 0..<f.len:
                for s in f[mu]:
                  var temp {.noinit.}: typeof(f[0][0])
                  temp := f[mu][s]*u[mu][s].adj
                  projectTAH(f[mu][s], temp)
          of Stout: discard
          of NoSmearing: self.setBoundaryConditions(u) # Undoes setting matter bc

  #[Optional force-gradient - not implemented]#
        
  # Update momentum
  if updateMomentum:
    threads:
      for mu in 0..<p.len:
        for s in p[mu]: p[mu][s] -= f[mu][s]

if isMainModule:
  qexInit()

  var 
    lat = intSeqParam("lat", @[4, 4, 4, 4])
    lo = lat.newLayout(@[1, 1, 1, 1])

    gaugeParams = %* {
      "action": "Wilson",
      "beta": 6.0,
      "steps": 10,
      "integrator": "MN2",
      "monte-carlo-algorithm": "hamiltonian-monte-carlo"
    }
    matterActionParams = %* {
      "smearing": "nHYP",
      "smearing-coefficients": @[0.4, 0.5, 0.5],
      "boundary-conditions": "aaaa"
    }
    fermionParams = %* {
      "mass": 0.0,
      "steps": 5,
      "integrator": "MN2",
      "monte-carlo-algorithm": "hamiltonian-monte-carlo"
    }
    bosonParams = %* {
      "mass": 0.1,
      "steps": 5,
      "integrator": "MN2",
      "monte-carlo-algorithm": "hamiltonian-monte-carlo"
    }

    rng = lo.newParallelRNG("MILC", 987654321)
    u = lo.newGauge()
    f = lo.newGauge()
    p = lo.newGauge()

  var gauge = lo.newGaugeAction(gaugeParams)
  var matter = lo.newMatterAction(matterActionParams)

  matter.addStaggeredFermion(lo, fermionParams)
  matter.addStaggeredBoson(lo, bosonParams)
  matter.addStaggeredHasenbuschFermion(lo, @[fermionParams, bosonParams, bosonParams])

  unit(u)
  unit(p)

  echo gauge.action(u)
  echo gauge.action()
  gauge.fields["gauge"].dtau = @[1.0]
  gauge.updateMomentum(u, f, p)

  echo matter.action(u, parallelRNG = some(rng))
  echo matter.action(u)
  for key, _ in matter.fields: matter.fields[key].dtau = @[1.0]
  matter.updateMomentum(u, f, p)

  qexFinalize()

proc updateGauge[A](self: LatticeFieldTheory; a: typedesc[A]) =
  case self.nested:
    of true: 
      var 
        subActions = newSeq[A]()
        nFields = 0
      for action in self.actions:
        for sAction in action.subActions:
          if not sAction.solo:
            if sAction.includeInStep: subActions.add sAction
      subActions.trajectory(0.5*self.tdtau)
      for sAction in subActions: sAction.includeInStep = false
    of false: self.u[].updateGauge(self.p, self.tdtau)

proc updateGauge*(self: LatticeFieldTheory) =
  self.updateGauge(self.AA)

proc updateMomentum*(self: LatticeFieldTheory; action: LatticeAction) =
  case action.action
    of PureGauge,GaugeMatter:
      let (updateGauge, gdtau) = self.f.getGaugeForce(action)
      case updateGauge:
        of true: self.updateMomentum(dtau = gdtau)
        of false: discard
    of PureMatter: discard
  
  case action.action:
    of PureGauge: discard
    of GaugeMatter,PureMatter:
      action.smear.rephased = false
      action.smear.smeared = false
      self.f.zero

      var updateStaggeredMatter = false
      for sAction in action.subActions:
        if self.f.getStaggeredMatterForce(sAction): 
          updateStaggeredMatter = true
      if updateStaggeredMatter: self.f.rephaseForce

      var updateWilsonMatter = false

      if (updateStaggeredMatter) or (updateWilsonMatter):
        self.f.setMatterBoundaryConditions(action.smear.bc)
        action.smear.sf(self.f, self.f)
        self.f.projectTrclssAntiHrm(self.u[])
        self.updateMomentum

var 
      step = 0
      pTime = 0.0
      cTime = 0.0
      pdtau = 0.0

    while true:
      var 
        minabsdt = Large64
        nNext = 0
        nRunning = 0
        updateSpace = false

      # Get minimum absolute step size & set V update step accordingly
      for action in self.actions:
        let absdts = action.getStep(step)
        for absdt in absdts:
          if absdt <= minabsdt: minabsdt = absdt
      for action in self.actions: cTime = action.setStep(step,minabsdt)
      self.tdtau = cTime - pTime + pdtau
      pTime = cTime

      # Update step
      for action in self.actions:
        nNext += action.getNumNext(step)
      if nNext == self.nFields: step += 1

      # Termination criterion
      for action in self.actions:
        nRunning += action.getNumRunning

      # Check if there is a "space" update
      for action in self.actions:
        for sAction in action.subActions:
          if abs(sAction.vdtau) > Small64: updateSpace = true

      # Either do "time/space" update or add current dtau to next dtau
      # Ensures that time update is not being done more than is necessary
      if (updateSpace) or (nRunning == 0):
        self.updateGauge
        for action in self.actions: self.updateMomentum(action)
        pdtau = 0.0
      else: pdtau = self.tdtau

      # Finish up
      if nRunning == 0: 
        for action in self.actions:
          for sAction in action.subActions:
            sAction.running = true
            sAction.iStep = 0
            sAction.timeStep = 0
            sAction.spaceStep = 0
            sAction.cumSum = 0.0
            sAction.absCumSum = 0.0
        break