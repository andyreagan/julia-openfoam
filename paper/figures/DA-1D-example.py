#!/usr/bin/python

import numpy as np

if __name__ == '__main__':
  from sys import argv

  ## essential plotting tool
  import matplotlib.pyplot as plt
  xkcd = False
  
  if argv[1] == 'xkcd':
    fig = plt.figure()
    plt.xkcd()
    xkcd = True
  else:
    ## use latex
    from matplotlib import rc
    rc('text', usetex=True)
    fig = plt.figure()

  

  xvec = np.linspace(-4,8,100)
  ## print xvec
  sigmaf = 1.5
  sigmao = 2.3
  uf = 0.0
  uo = 2.0
  ufdist = 1.0/(sigmaf*np.sqrt(2.0*np.pi))*np.exp(-((xvec-uf)**2)/(2*sigmaf**2))
  uodist = 1.0/(sigmao*np.sqrt(2.0*np.pi))*np.exp(-((xvec-uo)**2)/(2*sigmao**2))
  ## print ufdist

  ax1 = fig.add_axes([.2,.2,.7,.7])
  
  ## plot a line at0 (plot first, underneath timeseries)
  ax1.plot(xvec,ufdist,'b--')
  plt.ylim([0,.35])

  plt.yticks(fontsize=14)
  plt.xticks(fontsize=14)

  if xkcd:
    ax1.legend(['P(u(f))'],fontsize=16)
    plt.xlabel('u',fontsize=18)
    plt.ylabel('P(u)',fontsize=18)
    plt.savefig('DA-example01-labels-xkcd.pdf')
    plt.savefig('DA-example01-labels-xkcd.png')
  else:
    ax1.legend(['$P(u_f)$'],fontsize=16)
    plt.xlabel('$u$',fontsize=18)
    plt.ylabel('$P(u)$',fontsize=18)
    plt.savefig('DA-example01-labels.pdf')
    plt.savefig('DA-example01-labels.png')

  print 'figure 1 saved'

  ax1.plot(xvec,uodist,'g')
  if xkcd:
    ax1.legend(['P(u(f))','P(u(o))'])
    plt.savefig('DA-example02-labels-xkcd.pdf')
    plt.savefig('DA-example02-labels-xkcd.png')
  else:
    ax1.legend(['$P(u_f)$','$P(u_o)$'])
    plt.savefig('DA-example02-labels.pdf')
    plt.savefig('DA-example02-labels.png')

  print 'figure 2 saved'

  ## solve the 1D DA problem
  W = (sigmaf**2)/(sigmaf**2+sigmao**2)
  ## print W
  sigmaa = np.sqrt((1.0-W)*(sigmaf**2))
  ua = uf + W*(uo-uf)
  ## print ua

  ## plot the result
  uadist = 1.0/(sigmaa*np.sqrt(2*np.pi))*np.exp(-((xvec-ua)**2)/(2*sigmaa**2))
  ax1.plot(xvec,uadist,'k')

  if xkcd:
    ax1.legend(['P(u(f))','P(u(o))','P(u(a))'])
    plt.savefig('DA-example03-labels-xkcd.pdf')
    plt.savefig('DA-example03-labels-xkcd.png')
  else:
    ax1.legend(['$P(u_f)$','$P(u_o)$','$P(u_a)$'])
    plt.savefig('DA-example03-labels.pdf')
    plt.savefig('DA-example03-labels.png')

  print 'figure 3 saved'
  print 'success'



