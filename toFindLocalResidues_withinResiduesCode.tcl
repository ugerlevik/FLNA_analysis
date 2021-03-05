proc uniqueList {list} {
  set new {}
  foreach item $list {
    if {[lsearch $new $item] < 0} {
      lappend new $item
    }
  }
  return $new
}

set a [atomselect top "within 10 of resid 39 frame 0"]
uniqueList [$a get {resname resid chain}]
