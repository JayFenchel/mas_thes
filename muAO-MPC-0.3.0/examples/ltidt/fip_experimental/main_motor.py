r"""
WARNING: the fixed-point feature is experimental. It does not
check for overflow.  The computed resulst may be unprecise or
simply wrong.   

Example use of the experimental code generation of fixed-point 
version.  This should work on any target. 

The file main_motor.c included in this directory exemplifies
the use of the code. It automatically uses either fixed-point
or single precision float, according to the generated code.
"""
import muaompc
fixed_point = True 
mpc = muaompc.ltidt.setup_mpc_problem('sys_legomotor')
if fixed_point:
    mpc.generate_c_files(numeric='fip', fracbits=22)
else:
    mpc.generate_c_files(numeric='float32')
