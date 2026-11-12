#ctypes demo, this time from a C coded library
from ctypes import *
libc=CDLL.LoadLibrary("ident05_relsqlib.so.0")
create_tfdisc_c = libc.create_struct_tfdiscc

class IDENT05_TFDISC_C(Structure):
    _fields_ = [("nc", c_int),
                ("nd", c_int),
                ("size", c_int),
                ("Ts", c_double),
                ("coeffs", (c_double*4)),
                ("yinit", (c_double*2))]

create_tf_disc.argstype = [c_int, c_int, c_int, c_double, c_double*4, c_double*2] 
create_tf_disc.restype = IDENT05_TFDISC_C

FourPointsArrayType = c_double*4
TwoPointsArrayType = c_double*2

coeffs_t = FourPointsArrayType(0.025, 0.0, 0.975, 0.0)
yinit_t = TwoPointsArrayType(1.0, 0.5)
tfdisc_sc = create_tf_disc(1,1,4,0.1,coeffs_t, yinit_t)

# advance_c = libc.advanceintime
# advance_c.argstype=[POINTER(IDENT05_TFDISC_C), c_int)
# advance_c.restype = c_double
# yout = advance_c(tfdisc_sc, 1)