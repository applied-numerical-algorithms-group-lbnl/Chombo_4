#!/usr/bin/csh
./main.exe _expected/_case1/case1.inputs; mv pout.0 case1.pout.0
./main.exe _expected/_case2/case2.inputs; mv pout.0 case2.pout.0

diff case1.pout.0 _expected/_case1/mpi.1proc.2d.debug.spencer.pout.0 > diff_case1.out
diff case2.pout.0 _expected/_case2/mpi.1proc.2d.debug.spencer.pout.0 > diff_case2.out
