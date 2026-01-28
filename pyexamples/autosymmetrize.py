#!/usr/bin/python3 -i
import sys
sys.path = ['./']+sys.path
print(sys.path)
import qpp
from math import pi
#construxt Oh group
G=qpp.array_point_group()
G.generate(qpp.pq.RotMtrx(qpp.vector3(1,0,0),pi/2))
G.generate(qpp.pq.RotMtrx(qpp.vector3(0,1,0),pi/2))
G.generate(qpp.pq.RotMtrx(qpp.vector3(0,0,1),pi/2))
G.generate(qpp.pq.Sigma(qpp.vector3(0,0,1)))
len(G)
#use Oh group as generalied cell
C=qpp.pq.point_group_d([])
qpp.pq.generator_form(C,G)
C.end()
g=qpp.geometry(3)
g.set_cell(C)
g.typetable
g.set_cell(C)
g.build_types()
g.typetable.auto_update = True
g.typetable.auto_symmetrize = True
g.typetable.default_symmetrize_radius
g.typetable.default_symmetrize_radius=1.0
g.typetable.default_symmetrize_radius
g.add('U',.1,.1,.1)
g.add('O',1.8,.1,.1)
g.add('H',1.8,1.9,2.0)
