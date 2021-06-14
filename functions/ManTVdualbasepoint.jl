using Manopt

function ManTVdualbasepoint(m::PowPoint)
  d = size(m)
  dims = length(d)
  if dims > 1
    d2 = fill(1,dims+1)
    d2[dims+1] = dims
  else
    d2 = 1
  end
  n = PowPoint(repeat(getValue(m),inner=d2))

  return n

end
