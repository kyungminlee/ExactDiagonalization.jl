function make_bitmask(msb ::Integer;
                      lsb ::Integer=0,
                      dtype ::DataType=UInt)
  mask = dtype(0x1) << msb - dtype(0x1)
  submask = dtype(0x1) << lsb - dtype(0x1)
  return mask ^ submask
end