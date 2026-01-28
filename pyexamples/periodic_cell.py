import qpp
'''
This example showshow to replicate lattice nodes using periodic_cell and qpp.index
'''
qpp.REAL = 'float'
a = qpp.vector3(0,2.73,2.73)
b = qpp.vector3(2.73,0,2.73)
c = qpp.vector3(2.73,2.73,0)

cl = qpp.periodic_cell(a,b,c)

for I in qpp.index_range([-2,-2,-2],[2,2,2]):
    print(cl.transform(qpp.vector3(0),I))
