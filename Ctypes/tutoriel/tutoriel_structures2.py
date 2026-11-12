from ctypes import *

cdll.LoadLibrary("libc.so.6")

libc = CDLL("libc.so.6")

FourDouble = c_double * 4

class IDENT05_TESTDISC(Structure):
    _fields_ = [("Ts", c_int),
                ("size", c_int),
                ("maxsize", c_int),
                ("nc", c_int),
                ("nd", c_int),
                ("coeffs", POINTER(c_double))]

dcoeffs = FourDouble(1, -0.975, 0.025, 0.)
ex_test_disc = IDENT05_TESTDISC(100, 2, 2, 0.1, dcoeffs)

##ident05_testdisc = libc.ident05_testdisc
#printf.argtypes = [c_int, c_int, c_int, c_double, c_double, POINTER(c_double)]
 
