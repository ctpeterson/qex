#import base
#import stdUtils
#import comms
#import complexConcept
#import matrixConcept
#import metaUtils
import base/wrapperTypes
import macros
import typetraits

template makeDeclare(s:untyped):untyped {.dirty.} =
  template `declare s`*(t:typedesc):untyped {.dirty.} =
    template `declared s`*(y:t):untyped {.dirty.} = true
  template `is s`*(x:typed):untyped {.dirty.} =
    when compiles(`declared s`(x)):
      `declared s`(x)
    else:
      false
#makeDeclare(Scalar)
makeDeclare(Matrix)
makeDeclare(Vector)
makeDeclare(Real)
makeDeclare(Imag)
makeDeclare(Complex)

template forwardFunc*(t: typedesc, f: untyped) {.dirty.} =
  template f*(x: t): untyped =
    mixin f
    f(x[])


type
  Deref[T] = distinct T
template `[]`*[T](x: Deref[T]): untyped = (T(x))[]
template toDerefPtr*[T](x: T): untyped =
  when T is ptr or T is ref:
    x
  elif compiles(unsafeAddr x):
    Deref[ptr T](unsafeAddr x)
  else:
    x
template `[]`*(x: Deref, i: typed): untyped = x[][i]
template `[]`*(x: Deref, i,j: typed): untyped = x[][i,j]
forwardFunc(Deref, len)
forwardFunc(Deref, nrows)
forwardFunc(Deref, ncols)

#type
#  AsVar*[T] = object
#    v*: T
#template asVar*[T](x: T): untyped =
#  AsVar[type(T)](v: x)
#template asVar*(x: typed): untyped =
#  var tAsVar = x
#  tAsVar
#makeDeref(AsVar, x.T)
makeWrapperType(AsVar)
template len*(x: AsVar): untyped = x[].len
template re*(x: AsVar): untyped = x[].re
template im*(x: AsVar): untyped = x[].im
template `[]`*(x: AsVar; i: SomeInteger): untyped = x[][i]
template `[]`*(x: AsVar; i,j: SomeInteger): untyped = x[][i,j]
template assign*(r: AsVar, x: untyped) =
  var t = r[]
  t := x
template assign*(r: not AsVar, x: AsVar) =
  r := x[]
template `:=`*(r: AsVar, x: untyped) =
  var t = r[]
  t := x
template `*=`*(r: AsVar, x: untyped) =
  var t = r[]
  t *= x
template imul*(r: AsVar, x: untyped) =
  mixin imul
  var t = r[]
  imul(t, x)
forwardFunc(AsVar, nrows)
forwardFunc(AsVar, ncols)
forwardFunc(AsVar, numberType)
forwardFunc(AsVar, nVectors)
forwardFunc(AsVar, simdType)
forwardFunc(AsVar, simdLength)
forwardFunc(AsVar, norm2)
template norm2*(r: var auto, x: AsVar): untyped =
  mixin norm2
  norm2(r, x[])
template redot*(x: AsVar, y: AsVar): untyped =
  mixin redot
  redot(x[], y[])

makeWrapperType(Scalar)
template `[]`*(x: Scalar; i: typed): untyped = x[]
template `[]`*(x: Scalar; i,j: typed): untyped = x[]
#[
type
  Scalar*[T] = object
    v*: T
template asScalar*(x: typed): untyped =
  let xx = x
  AsScalar[type(xx)](v: xx)
makeDeref(AsScalar, x.T)
template `[]`*(x: AsScalar; i: SomeInteger): untyped = x[][i]
template `[]`*(x: AsScalar; i,j: SomeInteger): untyped = x[][i,j]
forwardFunc(AsScalar, nrows)
forwardFunc(AsScalar, ncols)
forwardFunc(AsScalar, numberType)
forwardFunc(AsScalar, nVectors)
forwardFunc(AsScalar, simdType)
forwardFunc(AsScalar, simdLength)
]#

#type
#  Adjointed*[T] = object
#    v*: T
#template adjointed*(x: typed): untyped =
  #static: echo "adjointed"
  #dumpTree: x
#  let x_adjointed = x
#  Adjointed[type(x_adjointed)](v: x_adjointed)
makeWrapperType(Adjointed)
template adj*(x: typed): untyped =
  mixin adj, isWrapper, asWrapper
  bind asAdjointed
  when isWrapper(x):
    #static: echo "adj typed wrapper"
    #dumpTree: x
    asWrapper(x, adj(x[]))
  elif compiles(adjImpl(x)):
    adjImpl(x)
  else:
    #static: echo "adj typed not wrapper"
    #dumpTree: x
    #(Masked[type(x)])(maskedObj(x,msk))
    asAdjointed(x)
#template `[]`*[T](x:Adjointed[T]):untyped = cast[T](x)
#makeDeref(Adjointed, x.T)
template `[]`*(x:Adjointed; i:SomeInteger):untyped = x[][i].adj
template `[]`*(x:Adjointed; i,j:SomeInteger):untyped = x[][j,i].adj
template len*(x:Adjointed):untyped = x[].len
template nrows*(x:Adjointed):untyped = x[].ncols
template ncols*(x:Adjointed):untyped = x[].nrows
template declaredVector*(x:Adjointed):untyped = isVector(x[])
template declaredMatrix*(x:Adjointed):untyped = isMatrix(x[])
template re*(x:Adjointed):untyped = x[].re
template im*(x:Adjointed):untyped = -(x[].im)
template simdType*(x: Adjointed): untyped = simdType(x[])
#template mvLevel*(x:Adjointed):untyped =
#  mixin mvLevel
#  mvLevel(x[])


type
  #ToSingle*{.borrow: `.`.}[T] = distinct T
  #ToSingle*[T] = distinct T
  ToSingle*[T] = object
    v*:T
template toSingleDefault*(xx: typed): untyped =
  lets(x,xx):
    when compiles(addr(x)):
    #when compiles(unsafeAddr(x)):
      #cast[ptr ToSingle[type(x)]](addr(x))[]
      cast[ptr ToSingle[type(x)]](unsafeAddr(x))[]
      #cast[ToSingle[type(x)]](x)
    else:
      #(ToSingle[type(x)])(x)
      cast[ToSingle[type(x)]](x)
      #cast[ToSingle[type(x)]]((var t=x; t.addr))
template toSingle*(xx: not typedesc):untyped =
  mixin isVector, isMatrix
  lets(x,xx):
    when compiles(toSingleImpl(x)):
      toSingleImpl(x)
    elif isComplex(x):
      asComplex(toSingle(x[]))
    elif isVector(x):
      asVector(toSingle(x[]))
    elif isMatrix(x):
      asMatrix(toSingle(x[]))
    elif x is SomeNumber:
      float32(x)
    else:
      toSingleDefault(x)
#template `[]`*[T](x:ToSingle[T]):untyped = cast[T](x)
makeDeref(ToSingle, x.T)
template `[]`*(x:ToSingle; i:SomeInteger):untyped = x[][i].toSingle
template `[]`*(x:ToSingle; i,j:SomeInteger):untyped = x[][j,i].toSingle
template len*(x:ToSingle):untyped = x[].len
template nrows*(x:ToSingle):untyped = x[].ncols
template ncols*(x:ToSingle):untyped = x[].nrows
template declaredVector*(x:ToSingle):untyped = isVector(x[])
template declaredMatrix*(x:ToSingle):untyped = isMatrix(x[])
template re*(x: ToSingle): untyped = toSingle(x[].re)
template im*(x: ToSingle): untyped = toSingle(x[].im)
template simdType*(x: ToSingle): untyped = simdType(x[])

type
  #ToDouble*{.borrow: `.`.}[T] = distinct T
  #ToDouble*[T] = distinct T
  ToDouble*[T] = object
    v*:T
#template toDoubleDefault*(xx: typed): untyped =
#  lets(x,xx):
#    when compiles(addr(x)):
#    #when compiles(unsafeAddr(x)):
#      #cast[ptr ToDouble[type(x)]](addr(x))[]
#      cast[ptr ToDouble[type(x)]](unsafeAddr(x))[]
#      #cast[ToDouble[type(x)]](x)
#    else:
#      #(ToDouble[type(x)])(x)
#      cast[ToDouble[type(x)]](x)
#      #cast[ToDouble[type(x)]]((var t=x; t.addr))
#template toDouble*(xx: typed): untyped =
#  mixin isVector, isMatrix, isComplex, toDoubleImpl
#  lets(x,xx):
#    when compiles(toDoubleImpl(x)):
#      toDoubleImpl(x)
#    elif isComplex(x):
#      asComplex(toDoubleDefault(x[]))
#    elif isVector(x):
#      asVector(toDouble(x[]))
#    elif isMatrix(x):
#      asMatrix(toDouble(x[]))
#    elif x is SomeNumber:
#      float64(x)
#    else:
#      toDoubleDefault(x)
template toDoubleX*[T](x: T): untyped =
  ToDouble[T](v: x)
template toDouble*(x: typed): untyped =
  mixin toDouble, toDoubleImpl, isWrapper, asWrapper
  when isWrapper(x):
    #static: echo "toDouble typed wrapper"
    #dumpTree: x
    asWrapper(x, toDouble(x[]))
  else:
    #static: echo "toDouble typed not wrapper"
    #dumpTree: x
    #(Masked[type(x)])(maskedObj(x,msk))
    toDoubleImpl(x)
#template `[]`*[T](x:ToDouble[T]):untyped = cast[T](x)
makeDeref(ToDouble, x.T)
#template `[]`*(x:ToDouble):untyped = x.v
#template `[]=`*(x:t; y:auto):untyped =
#   x.v = y
template `[]`*(x:ToDouble; i:SomeInteger):untyped = x[][i].toDouble
template `[]`*(x:ToDouble; i,j:SomeInteger):untyped = x[][j,i].toDouble
template len*(x:ToDouble):untyped = x[].len
template nrows*(x:ToDouble):untyped = x[].nrows
template ncols*(x:ToDouble):untyped = x[].ncols
template declaredVector*(x:ToDouble):untyped = isVector(x[])
template declaredMatrix*(x:ToDouble):untyped = isMatrix(x[])
template re*(x:ToDouble):untyped = toDouble(x[].re)
template im*(x:ToDouble):untyped = toDouble(x[].im)
template simdType*(x: ToDouble): untyped = simdType(x[])
macro dump2(x: typed): auto =
  result = newEmptyNode()
  echo x.treerepr
template numberType*(x: ToDouble): untyped =
  dump2: x
  numberType(x[])


template sameWrapper*(x: typedesc, y: typedesc): untyped =
  type(stripGenericParams(x)) is type(stripGenericParams(y))

#  FIXME: should get a pointer
template getAlias*(x: typed): untyped =
  when compiles(unsafeAddr(x)):
    unsafeAddr(x)
  else:
    var tGetAlias = x
    addr(tGetAlias)

# index I should be wrapped: Simd[int] instead of int
type
  Indexed*[T,I] = object
    indexedPtr*: T
    indexedIdx*: I

template indexedX[T,I](x: T, i: I): untyped =
  Indexed[T,I](indexedPtr: x, indexedIdx: i)
template indexed*[T,I](x: T, i: I): untyped =
  when isWrapper(T):
    #static: echo "indexed isWrapper"
    when sameWrapper(type T, type I):
      #static: echo "indexed sameWrapper"
      x[i[]]
    else:
      #static: echo "indexed not sameWrapper"
      #var tIndexed = asWrapper(T, indexedX(getAlias x[], i))
      var tIndexed = asWrapper(T, indexed(x[], i))
      tIndexed
  else:
    #static: echo "indexed not isWrapper"
    var tIndexed = indexedX(getAlias x, i)
    tIndexed
template mindexed*[T,I](x: T, i: I): untyped =
  when isWrapper(T):
    #static: echo "mindexed isWrapper"
    when sameWrapper(type T, type I):
      #static: echo "mindexed sameWrapper"
      x[i[]]
    else:
      #static: echo "mindexed not sameWrapper"
      #var tIndexed = asWrapper(T, indexedX(getAlias x[], i))
      var tIndexed = asWrapper(T, indexed(x[], i))
      tIndexed
  else:
    #static: echo "mindexed not isWrapper"
    var tIndexed = indexedX(getAlias x, i)
    tIndexed
template obj(x:Indexed):untyped = x.indexedPtr[]
template idx(x:Indexed):untyped = x.indexedIdx
template `!`(x:Indexed):untyped = obj(x)[idx(x)]
template `!=`(x:Indexed,y:typed):untyped = obj(x)[idx(x)] = y
template optIndexed*[T,I](x: T, i: I): untyped =
  x[i]
  #when isWrapper(T):
  #  #static: echo "optindexed isWrapper"
  #  when sameWrapper(type T, type I):
  #    #static: echo "optindexed sameWrapper"
  #    x[i[]]
  #  else:
  #    #static: echo "optindexed not sameWrapper"
  #    #var tIndexed = asWrapper(T, indexedX(getAlias x[], i))
  #    var tIndexed = asWrapper(T, indexed(x[], i))
  #    tIndexed
  #else:
  #  #static: echo "optindexed not isWrapper"
  #  x[i]
template `[]`*(x:Indexed):untyped =
  let tIndexedBracket0 = x
  obj(tIndexedBracket0)[idx(tIndexedBracket0)]
template `[]`*(x:Indexed; i:SomeInteger):untyped =
  let tIndexedBracket1 = x
  optIndexed(obj(tIndexedBracket1)[i], idx(tIndexedBracket1))
template `[]`*(x:Indexed; i,j:SomeInteger):untyped =
  let tIndexedBracket2 = x
  optIndexed(obj(tIndexedBracket2)[i,j], idx(tIndexedBracket2))
template `[]=`*(x:Indexed; i:SomeInteger; y: typed) =
  let tIndexedBracket1Eq = x
  optIndexed(obj(tIndexedBracket1Eq)[i], idx(tIndexedBracket1Eq)) := y
template `[]=`*(x:Indexed; i,j:SomeInteger; y: typed) =
  let tIndexedBracket2Eq = x
  #optIndexed(obj(tIndexedBracket2Eq)[i,j], idx(tIndexedBracket2Eq)) := y
  obj(tIndexedBracket2Eq)[i,j][idx(tIndexedBracket2Eq)] = y
template `:=`*(x:Indexed, y: typed) =
  #echoRepr: x
  #echoRepr: y
  let tIndexedColonEq = x
  #when isWrapper(type(obj(tIndexedColonEq))):
  #  mixin `:=`
  #  var tIndexedColonEqWrap = mindexed(obj(tIndexedColonEq),idx(tIndexedColonEq))
  #  tIndexedColonEqWrap := y
  #else:
  obj(tIndexedColonEq)[idx(tIndexedColonEq)] = y
template assign*(x: SomeNumber, y: Indexed) =
  x := y[]
template `+=`*(x:Indexed, y: typed) =
  let tIndexedPlusEq = x
  tIndexedPlusEq != !tIndexedPlusEq + y
template `-=`*(x:Indexed, y: typed) =
  let tIndexedMinusEq = x
  tIndexedMinusEq != !tIndexedMinusEq - y
template `*=`*(x:Indexed, y: typed) =
  let tIndexedStarEq = x
  tIndexedStarEq != !tIndexedStarEq * y

template len*(x:Indexed):untyped = obj(x).len
template nrows*(x:Indexed):untyped = obj(x).nrows
template ncols*(x:Indexed):untyped = obj(x).ncols
#template declaredVector*(x:Indexed):untyped = isVector(x.obj)
#template declaredMatrix*(x:Indexed):untyped = isMatrix(x.obj)
template re*(x:Indexed):untyped =
  let tIndexedRe = x
  obj(tIndexedRe).re[idx(tIndexedRe)]
template im*(x:Indexed):untyped =
  let tIndexedIm = x
  obj(tIndexedIm).im[idx(tIndexedIm)]
template `re=`*(x: Indexed, y: typed): untyped =
  let tIndexedReEq = x
  obj(tIndexedReEq).re[idx(tIndexedReEq)] = y
template `im=`*(x: Indexed, y: typed): untyped =
  let tIndexedImEq = x
  obj(tIndexedImEq).im[idx(tIndexedImEq)] = y
template simdType*(x: Indexed): untyped = simdType(obj(x))
template numberType*(x: Indexed): untyped = numberType(obj(x))
template `*`*(x: Indexed, y: typed): untyped =
  let tIndexedMul = x
  obj(tIndexedMul)[idx(tIndexedMul)] * y

# FIXME?: Masked, VarMasked
type
  MaskedObj*[T] = object
    pobj*: ptr T
    mask*: int
  Masked*[T] = MaskedObj[T]
  #Masked2*[T] = Masked[T]
#template pobj*(x:Masked):untyped = ((MaskedObj[x.T])(x)).pobj
#template mask*(x:Masked):untyped = ((MaskedObj[x.T])(x)).mask
#template `pobj=`*(x:Masked;y:untyped):untyped = ((MaskedObj[x.T])(x)).pobj = y
#template `mask=`*(x:Masked;y:untyped):untyped = ((MaskedObj[x.T])(x)).mask = y
template maskedObj*[T](x: T, msk: int): untyped =
  MaskedObj[type(T)](pobj: addr(x), mask: msk)
  #MaskedObj[type(T)](pobj: addr(x).regenSym, mask: msk.regenSym)
template maskedX*(x: typed, msk: int): untyped =
  bind maskedObj
  mixin isWrapper
  when isWrapper(x):
    asWrapper(x, maskedX(x[], msk))
  #elif x is SomeNumber:
  #  x
  else:
    maskedObj(x, msk)
template masked*(x: typed, msk: int): untyped =
  when x is SomeNumber:
    x
  else:
    bind maskedX
    var tMasked = maskedX(x, msk)
    tMasked
template varMasked*(x: typed, msk: int): untyped =
  when x is SomeNumber:
    x
  else:
    bind maskedX
    var tVarMasked = maskedX(x, msk)
    tVarMasked
template `[]`*(m: Masked): untyped = m.pobj[]
macro `[]`*(m: Masked{nkObjConstr}): untyped = m[1][1]
template getMask*(m: Masked): untyped = m.mask
template index1U*(m: Masked; i: typed): untyped =
  masked(m[][i], m.getMask)
  #let tMaskedIndex1 = m
  #masked(tMaskedIndex1[][i], tMaskedIndex1.mask)
template `[]`*(m: Masked; i: typed): untyped =
  flattenCallArgs(index1U, m, i)
template `[]`*(m: Masked; i,j: typed): untyped =
  let tMaskedIndex2 = m
  masked(tMaskedIndex2[][i,j], tMaskedIndex2.mask)
#template isRIC*(x:int):untyped = true
#template isRIC*(m:Masked):untyped = isRIC(m.pobj[])
#template isComplex*(m:Masked):untyped = isComplex(m.pobj[])
template declaredComplex*(m:Masked):untyped =
  mixin declaredComplex
  declaredComplex(m.pobj[])
template isVector*(m:Masked):untyped =
  mixin isVector
  isVector(m.pobj[])
#template isMatrix*(m:Masked):untyped =
#  mixin isMatrix
  #echo "isMatrix"
  #echo isMatrix(m.pobj[])
#  isMatrix(m.pobj[])
#template mvLevel*(m:Masked):untyped =
#  mixin mvLevel
#  mvLevel(m.pobj[])
template numNumbers*(m:Masked):untyped = numNumbers(m[])
template numberType*(m:Masked):untyped = numberType(m[])
template len*(m:Masked):untyped =
  mixin len
  len(m.pobj[])
template nrows*(m:Masked):untyped =
  mixin nrows
  nrows(m.pobj[])
template ncols*(m:Masked):untyped =
  mixin ncols
  ncols(m.pobj[])
template re*(m: Masked): untyped =
  mixin re
  let tMaskedRe = m
  masked(tMaskedRe.pobj[].re, tMaskedRe.mask)
template im*(m: Masked): untyped =
  mixin im
  let tMaskedIm = m
  masked(tMaskedIm.pobj[].im, tMaskedIm.mask)
#template assign*(m: Masked; x: SomeNumber): untyped =
#  m.pobj[] = x
template `re=`*(m: Masked; x: auto): untyped =
  mixin re
  let tMaskedReEq = m
  assign(masked(tMaskedReEq.pobj[].re, tMaskedReEq.mask), x)
template `im=`*(m: Masked; x: auto): untyped =
  mixin im
  let tMaskedImEq = m
  assign(masked(tMaskedImEq.pobj[].im, tMaskedImEq.mask), x)
proc assign*[T,U](m:Masked[T], x:Masked[U]) =
  ## Only works for the same number of unmasked bits,
  ## and assign those from RHS to LHS in sequence.
  var
    i,j = 0
    b = m.mask
    c = x.mask
  while b != 0:
    if (b and 1) != 0:
      while c != 0:
        let p = (c and 1) != 0
        if p: m.pobj[][i] = x.pobj[][j]
        c = c shr 1
        j.inc
        if p: break
    b = b shr 1
    i.inc
proc `:=`*[T,U](m:Masked[T], x:Masked[U]) = assign(m,x)
proc assign*[T](m: Masked[T], x: SomeNumber) =
  var
    i = 0
    b = m.mask
  while b != 0:
    if (b and 1) != 0:
      m.pobj[][i] = x
      break
    b = b shr 1
    i.inc
proc `:=`*[T](m:Masked[T], x: SomeNumber) = assign(m,x)
template defopeq(oe,o:untyped) =
  proc oe*[T](m:Masked[T], x: SomeNumber) =
    var t: type(x)
    t := m
    m := o(t,x)
defopeq(`+=`, `+`)
defopeq(`-=`, `-`)
defopeq(`*=`, `*`)
defopeq(`/=`, `/`)
proc assign*[T](m: var SomeNumber, x:Masked[T]) =
  var
    i = 0
    b = x.mask
  while b != 0:
    if (b and 1) != 0:
      m = x.pobj[][i]
      break
    b = b shr 1
    i.inc
template `:=`*[T](m: SomeNumber, x: Masked[T]) = assign(m,x)
template `+=`*[T](m: SomeNumber, x: Masked[T]) =
  var t: type(m)
  t := x
  m += t
proc norm2*(m: Masked): float =
  #var r: type(toDouble(m.pobj[][0]))
  var
    i = 0
    b = m.mask
  while b != 0:
    if (b and 1) != 0:
      let t = m.pobj[][i]
      result += t*t
    b = b shr 1
    i.inc

#template eval*(x: AsComplex): untyped = asComplex(eval(x[]))
template eval*(x: ToDouble): untyped =
  #echoType: x
  mixin map
  #map(map(x[],toDouble),eval)
  template etd(y: untyped): untyped = eval(toDouble(y))
  map(x[],etd)
template eval*(x: SomeNumber): untyped = x
#template eval*(x: typed): untyped =
#  mixin isComplex
#  when isComplex(x):
#    asComplex(eval(x[]))
#  elif x is SomeNumber:
#    x
#  else:
#    map(map(x[],toDouble),eval)
