mpirun -np 16 code/bin/athena -r BDstream.final.rst time/tlim=45.0 problem/ecc=0.2 problem/change_setup=1 job/problem_id=BDrst_ecc >ecc.log&

#mpirun -np 16 code/bin/athena -r BDstream.final.rst time/tlim=25.0 problem/ecc=0.0 problem/change_setup=1 job/problem_id=BDrst_circ >circ.log&


