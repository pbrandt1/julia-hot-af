nheads = @parallel (+) for i=1:100000000
  int(randbool())
end

nheads
