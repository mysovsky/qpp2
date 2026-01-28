#!/usr/bin/python3 -i
import sys
sys.path = ['./']+sys.path
print(sys.path)
import qpp
from math import pi
Td = qpp.pq.array_point_group_d()
Td.generate(qpp.pq.RotMtrx(qpp.vector3(1,1,1),2*pi/3))
Td.generate(qpp.pq.RotMtrx(qpp.vector3(-1,-1,1),2*pi/3))
Td.generate(qpp.pq.Sigma(qpp.vector3(1,1,0)))
a=2.73
v1=qpp.vector3(0,a,a)
v2=qpp.vector3(a,0,a)
v3=qpp.vector3(a,a,0)
bx = qpp.periodic_cell(v1,v2,v3)
