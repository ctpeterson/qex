type
  MCStream* = object
    name*: string
    output*: string

proc newMCStream*(name: string; start: bool = false): MCStream = 
  result = MCStream(name: name)
  echo "<begin: " & result.name & ">"

proc add*(self: var MCStream; text: string) =
  echo self.output & text & "\n"

proc finishStream*(self: var MCStream) =
  echo self.output & "<end: " & self.name & ">"