#!/usr/bin/python3 -i
import qpp
import random
pq.globals.ncores=6
N=2000
g=qpp.geometry(3)
g.cell[0]=qpp.vector3(1,0,0)
g.cell[1]=qpp.vector3(0,1,0)
g.cell[2]=qpp.vector3(0,0,1)
for i in range(N):
    g.add('X', random.random(),  random.random(),  random.random())
bt=qpp.bonding_table()
bt.default=.1
g.build_types()
nt=qpp.ngbr_table(g,bt)
nt.build()
