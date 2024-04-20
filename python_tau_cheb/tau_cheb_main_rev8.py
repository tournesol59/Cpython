import math as math
import lhlibraylin as lhlu # C extension package

#recopy of the MAIN part program present in tau_cheb_methodfine_rev8.py
#common part: control parameters
taun=0.09
M_order=7;

#third part: Fourier cospi
decump=calc_Fourier_cospi(4)  # bit strange but ok: cos2,cos3,cos6,cos12
c0=decump[0]
c1=decump[1]
d1=decump[2]
c32=decump[3]
d32=decump[4]
c3=decump[5]
d3=decump[6]
c6=decump[7]
d6=decump[8]
print("Fourier coeffs: "+str(c0)+" : "+str(c1)+" : "+str(d1)+" : "+str(c32)+" : "+str(d32)+" : "+str(c3)+" : "+str(d3)+" : "+str(c6)+" : "+str(d6))
# the following comes from the inversion x^i <-> Ti(x) of cos(3x), sin(3x) in Matlab
Ucos5=[   -0.3008,  0,  -1.0371, 0,  0.2320, 0]
Usin5=[   0, 0.8906, 0, -0.4922,  0, 0.1266]

#below: sign (-1) because the change of variable gives cos(2*pi*x)~= cos(pi*z+pi)= -cos(pi*z)
Uvar1 = [0.0 for i in range(M)]
Uvar1[0] =  (-(lamb/8.))*(-c32)*Ucos5[0]  # %(-(lamb/8.+pi^2/40))*(c3/2)*Ucos6(1:2:5) ;  /8 =1/(2^2)/2
Uvar1[0] = Uvar1[0] -(lamb/8)*(-c0) -lamb/8 #%Uvar1(1) -(lamb/8+pi^2/40)*(c0/2) -lamb/8 
Uvar1[1] = (-(lamb/8.))*(-d32)*Usin5[1] 
Uvar1[2] =  (-(lamb/8.))*(-c32)*Ucos5[2]
Uvar1[3] = (-(lamb/8.))*(-c32)*Usin5[3]
Uvar1[4] = taun
Uvar1[5] = 0.
Uvar1[6] = y_limit # % see true solution by Fourier series below to see that u(1) ~= -0.25
#export in screen
print("Uvar1")
print(str(Uvar1))

# ************** then SOLVING ************* in this version, try an own lu() limited to order of 7
Anum=main_subroutine_der(M)
print("LU decumposition of Anum")
print(str(Anum))
#call to our own C extension module
lhlu.lu_parser_wrapper("input",8);
#so far we need a file transfer to get the matrix result
file_res = open("result")
print(fileres.read())
file_res.close()

