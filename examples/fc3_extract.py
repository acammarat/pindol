import h5py
f = h5py.File("fc3.hdf5",'r')
g = f['fc3'][()]
f.close()
outf = "fc3.dat"
out = open(outf,"w")

print(" Extracting third-order force constants from file fc3.hdf5")

for i in range(g.shape[0]):
  for j in range(g.shape[0]):
    for k in range(g.shape[0]):
      for l in range(3):
        for m in range(3):
          for n in range(3):
            print("%25.22E" % (g[i][j][k][l][m][n]), file=out)
#            print("%25.22E" % (g[i][j][k][l][m][n]))
#            print(" %i %i %i %i %i %i" % (i,j,k,l,m,n))

print(" Output written in file %s." % outf)

out.close()

