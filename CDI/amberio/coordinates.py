#! /usr/bin/env python
################################################################################
#                                                                              
# FILE: coordinates.py - 
#
# DESCRIPTION:
#
# AUTHOR: Brian K. Radak (BKR) - <radakb@biomaps.rutgers.edu>
#
################################################################################
from math import sqrt,acos,pi
def Bond(crds,i,j):
    return BondFromVecs(crds[3*i:3*(i+1)],crds[3*j:3*(j+1)])
    
def BondAndGradients(crds,i,j):
    pass

def BondFromVecs(a,b):
    return VecMag(VecDiff(a,b))

def Angle(crds,i,j,k):
    return AngleFromVecs(crds[3*i:3*(i+1)],crds[3*j:3*(j+1)],crds[3*k:3*(k+1)])

def AngleAndGradients():
    pass

def AngleFromVecs(a,b,c):
    rBA = BondFromVecs(b,a)
    rBC = BondFromVecs(b,c)
    rAB = BondFromVecs(a,b)
    return acos( (rBA**2 + rBC**2 - rAB**2) / (2.*rBA*rBC) )

def Dihedral(crds,i,j,k,l):
    vi = crds[3*i:3*(i+1)]
    vj = crds[3*j:3*(j+1)]
    vk = crds[3*k:3*(k+1)]
    vl = crds[3*l:3*(l+1)]
    return DihedralFromVecs(vi,vj,vk,vl)

def DihedralAndGradients():
    pass

def DihedralFromVecs(a,b,c,d):
    vAB = VecDiff(a,b)
    vBC = VecDiff(b,c)
    vCD = VecDiff(c,d)
    nABC = VecCross(vAB,vBC)
    nBCD = VecCross(vBC,vCD)
    denom = VecMag(nABC)*VecMag(nBCD)
    dihedral = 0.
    if denom > 1.e-30:
        cosphi = -VecDot(nABC,nBCD)/denom
        if cosphi <= -1.:  dihedral = pi
        elif cosphi >= 1.: dihedral = 0.
        else:              dihedral = acos(cosphi)
        if VecDot(vBC,VecCross(nABC,nBCD)) > 0.: dihedral *= -1.
        dihedral = pi - dihedral
    return dihedral

def VecMag(u):
    return sqrt(VecDot(u,u))

def VecSum(u,v):
    return [ ui + vi for ui,vi in zip(u,v) ]

def VecDiff(u,v):
    return [ ui - vi for ui,vi in zip(u,v) ]

def VecDot(u,v):
    dot = 0.
    for ui,vi in zip(u,v): dot += ui*vi
    return dot

def VecCross(u,v):
    return [ u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0] ]

#
# TODO: Redo this to get gradients?
# 
# def Dihedral(Ra,Rb,Rc,Rd):
#     d = 0.
#     Rba = VecDiff(Rb,Ra)
#     Rcb = VecDiff(Rc,Rb)
#     Rdc = VecDiff(Rd,Rc)
#     #Rca = VecDiff(Rc,Ra) # only needed for gradient
#     #Rdb = VecDiff(Rd,Rb) # only needed for gradient
 
#     A = np.cross(Rba,Rcb)
#     B = np.cross(Rcb,Rdc)
#     n2A = (A**2).sum()
#     n2B = (B**2).sum()
#     den = np.sqrt(n2A*n2B)


#     print "den",den
#     if den > 1.e-30:
#         x = -(A[0]*B[0] + A[1]*B[1] + A[2]*B[2])/den;
#         print "x",x
#         if x <= -1.:
#             d = np.pi
#         elif x >= 1.:
#             d = 0.
#         else:
#             d = np.arccos( x )

#         BxA = np.cross(B,A)

#         BondProj = Rcb[0]*BxA[0] + Rcb[1]*BxA[1] + Rcb[2]*BxA[2]
#         if BondProj > 0.0: 
#             d = -d
#         d = np.pi - d

#         # ncb = math.sqrt(Rcb[0]*Rcb[0]+Rcb[1]*Rcb[1]+Rcb[2]*Rcb[2])

#         # sclA = 1. / ( n2A*ncb )
#         # A[0] *= sclA
#         # A[1] *= sclA
#         # A[2] *= sclA
#         # AxRcb = np.cross(A,Rcb)

        # sclB = 1. / ( n2B*ncb )
        # B[0] *= sclB
        # B[1] *= sclB
        # B[2] *= sclB
      
        # RcbxB = np.cross(Rcb,B)
        # ddda = np.cross(AxRcb,Rcb)
        # dddb = np.cross(Rca,AxRcb)
        # foo = np.cross(RcbxB,Rdc)
        # dddb = [ dddb[i] + foo[i] for i in range(3) ]
        # dddc = np.cross(AxRcb,Rba)
        # foo = np.cross(Rdb,RcbxB)
        # dddc = [ dddc[i] + foo[i] for i in range(3) ]
        # dddd = np.cross(RcbxB,Rcb)

#    return d
