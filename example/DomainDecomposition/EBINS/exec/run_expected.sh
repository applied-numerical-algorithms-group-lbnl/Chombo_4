#!/usr/bin/csh
./main.exe _expected/_case1/inflow_outflow_treb.inputs     |& tee case1.out
./main.exe _expected/_case2/inflow_outflow_eb.inputs       |& tee case2.out
./main.exe _expected/_case3/conservative_gradient.inputs   |& tee case3.out

diff case1.out _expected/_case1/2d.debug.serial.out > diff_case1.out
diff case2.out _expected/_case2/2d.debug.serial.out > diff_case2.out
diff case3.out _expected/_case3/2d.debug.serial.out > diff_case3.out
