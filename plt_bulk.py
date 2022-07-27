#!/usr/bin/env python3

import numpy as np

def bulk_disp(h,k,l,branch,param,S):

  DzA  = param[0]
  DzB  = param[1]
  J1   = param[2]
  J2   = param[3]
  A2   = param[4]
  J3   = param[5]
  Jp1  = param[6]
  Jp2A = param[7]
  Jp2B = param[8]
  Jp3A = param[9]
  Jp3B = param[10]

  U = np.exp( 2j*np.pi * h ) + np.exp( 2j*np.pi * k ) + np.exp (2j*np.pi * l ) 
  V = np.exp( 2j*np.pi * (h-k+l) ) + np.exp( 2j*np.pi * (-h+k+l) ) + np.exp( 2j*np.pi * (h+k-l) ) 
  W = np.exp( 2j*np.pi * (h+l) ) + np.exp( 2j*np.pi * (k+l) ) + np.exp( 2j*np.pi * (h+k) ) 
  X = np.exp( 2j*np.pi * (l-h) ) + np.exp( 2j*np.pi * (k-h) ) + np.exp( 2j*np.pi * (k-l) ) 

  DMnn = 2.* ( np.sin( 2*np.pi * (h-l) ) - np.sin( 2*np.pi * (k-l) ) - np.sin( 2*np.pi * (h-k) ))

  JAA = J2 *2.*np.real(X) + Jp2A*2.*np.real(U) + A2 * DMnn 
  JBB = J2 *2.*np.real(X) + Jp2A*2.*np.real(U) - A2 * DMnn

  IA = 3*J1 + 6*J2 + 3*J3 + 2*DzA + Jp1 + 6*Jp2A + 3*Jp2B + 3*Jp3A + 3*Jp3B
  IB = 3*J1 + 6*J2 + 3*J3 + 2*DzB + Jp1 + 6*Jp2A + 3*Jp2B + 3*Jp3A + 3*Jp3B

  JAB = J1*U + J3*V + Jp1 + Jp2B*W + Jp3A*X + Jp3B*np.conj(X)

  AA = IA - JAA
  BB = IB - JBB
  AB = JAB


  if branch == 0:
     res = -0.5*S*( (AA+BB) + np.sqrt( (AA+BB)**2. - 4.*(AA*BB - np.asb(AB)**2.) )) 
  else:
     res = -0.5*S*( (AA+BB) - np.sqrt( (AA+BB)**2. - 4.*(AA*BB - np.asb(AB)**2.) )) 
  return res



