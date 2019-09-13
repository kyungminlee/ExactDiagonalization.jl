function make_bitmask{U}(msb ::Integer;
                         lsb ::Integer=0) ::U where {U <: Unsigned}
    mask = U(0x1) << msb - U(0x1)
    submask = U(0x1) << lsb - U(0x1)
    return mask ^ submask
end