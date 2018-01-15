import ../comms/mpi
import posix
import times
import os
import strutils

type TimerObj = tuple[s:string,d:float,n:int]
var timer: float
var deltas = newSeq[ptr TimerObj](0)
template tic =
  timer = epochTime()
template toc(st: string) =
  block:
    var first{.global.} = true
    var to{.global.}: TimerObj
    if first:
      first = false
      deltas.add to.addr
      to.s = st
      to.d = 0.0
      to.n = 0
    let t = epochTime()
    to.d += t - timer
    timer = t
    inc to.n

discard MPI_init()
let comm = MPI_COMM_WORLD
var rank, size: cint
discard MPI_Comm_rank(comm, rank.addr)
discard MPI_Comm_size(comm, size.addr)
var fn = "testout.dat"
if paramCount()>=1:
  fn = paramStr(1)
var sinfo = true

proc testContig(n: int) =
  var wbuf = newSeq[cint](n)
  var rbuf = newSeq[cint](n)
  for i in 0..<n:
    wbuf[i] = cint(i+1)

  let wflags = O_CREAT or O_WRONLY
  let wmode = 0o666
  let wsize = wbuf.len*sizeof(type(wbuf[0]))
  let woffset = rank*wsize

  let rflags = O_RDONLY
  let rmode = 0o666
  let rsize = rbuf.len*sizeof(type(rbuf[0]))
  let roffset = rank*rsize

  tic()
  let wfd = open(fn, wflags, wmode)
  discard MPI_barrier(comm)
  toc("open write")
  let woff2 = lseek(wfd, woffset, SEEK_SET)
  toc("seek write")
  let wsize2 = write(wfd, wbuf[0].voidaddr, wsize)
  discard MPI_barrier(comm)
  toc("write")
  discard close(wfd)
  toc("write close")

  tic()
  let rfd = open(fn, rflags, rmode)
  discard MPI_barrier(comm)
  toc("open read")
  let roff2 = lseek(rfd, roffset, SEEK_SET)
  toc("seek read")
  let rsize2 = read(rfd, rbuf[0].voidaddr, rsize)
  discard MPI_barrier(comm)
  toc("read")
  discard close(rfd)
  toc("read close")

  var errors,errors2 = 0
  for i in 0..<n:
   if rbuf[i] != wbuf[i]:
     inc errors
     break
  discard MPI_Allreduce(errors.voidaddr, errors2.voidaddr,
                        1.cint, MPI_INT64_T, MPI_SUM, comm)
  if errors2>0:
    if rank==0:
      echo "errors: ", errors2 / size

  if rank==0: removeFile(fn)


proc disp =
  if rank==0:
    echo "total ranks: ", size
    echo "using file: ", fn
    echo "sizeof(Off): ", sizeof(Off)
    {.emit:"""printf("sizeof(off_t): %i\n", sizeof(off_t));""".}
disp()

var n = 1024
var nmax = 2*72*16*1024
if paramCount()>=2:
  nmax = parseInt(paramStr(2))*1024*1024
while n<=nmax:
  if rank==0:
    let b = formatSize(n*size*sizeof(cint))
    echo "total file size: ", b
  testContig(n)
  sinfo = false
  testContig(n)
  for i in 0..<deltas.len:
    deltas[i].d = 0.0
    deltas[i].n = 0
  testContig(n)
  testContig(n)
  for i in 0..<deltas.len:
    if rank==0:
      let sp = max(0, 23-deltas[i].s.len)
      let t = formatFloat(1e6*deltas[i].d/deltas[i].n.float, ffDecimal, 3)
      echo deltas[i].s, spaces(sp), ": ", align(t,15)
  n *= 2

discard MPI_Finalize()
