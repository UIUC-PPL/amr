#
import numpy as na
import os
import matplotlib
import matplotlib.pyplot as plt

class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value
class Plot:
	def __init__(self,grids):
		plt.figure(1,figsize=(8,8))
		for b in grids:
			nx = b["nx"]; ny = b["ny"]
			dx = b["dx"]; dy = b["dy"]
			plt.imshow(b["u"],extent=(b["xmin"]-0.49*dx,b["xmax"]+0.49*dx,
					b["ymin"]-0.49*dy,b["ymax"]+0.49*dy))
			plt.clim([1.,2.])
		plt.xlim([0.,1.])
		plt.ylim([0.,1.])
		plt.xlabel("X")
		plt.ylabel("Y")
		plt.savefig("out1.png")
		plt.show()
		return

class Load:
	def __init__(self,step,dataPath="./out"):
		scriptPath = os.getcwd()
		if dataPath==None:
			dataPath=scriptPath
		self._dataPath   = dataPath
                self._scriptPath = scriptPath
		blocks = self.loadData(step)
		grids  = self.constructGrids(blocks)
		self.grids = grids
		return 
	def constructGrids(self,blocks):
		grids = []
		for i,b in enumerate(blocks):
			grid = AutoVivification()
			xpos = b["xpos"]
			ypos = b["ypos"]
			value= b["u"]
			xmin = min(xpos)
			ymin = min(ypos)
			xmax = max(xpos)
			ymax = max(ypos)
			grid["nx"] = na.sqrt(len(xpos))
			grid["ny"] = na.sqrt(len(ypos))
			xsize = (xmax-xmin)/(grid["nx"]-1.)*grid["nx"]; dx = xsize/grid["nx"]
			ysize = (ymax-ymin)/(grid["ny"]-1.)*grid["ny"]; dy = ysize/grid["ny"]
			grid["xmin"] = xmin
			grid["ymin"] = ymin 
			grid["xmax"] = xmax
			grid["ymax"] = ymax
			grid["index"] = i
			nx = grid["nx"] 
			ny = grid["ny"] 
			u = na.zeros(nx*ny).reshape(nx,ny)

			n = 0
			x = xpos[0:-1:nx]
			y = ypos[0:ny]
			grid["dx"] = x[1]-x[0]
			grid["dy"] = y[1]-y[0]

			for i in range(int(nx)):
				for j in range(int(ny)):	
					u[ny-1-j][i] = value[n]
					n+=1
			grid["u"] = u
			grids.append(grid)
			del grid
		return grids
	def isDataFiles(self,s):
		if s.find("out_")==-1:
			return False
		else:
			if s.find("swp")==-1:
				return True
			else:
				return False
	def loadData(self,step):
		if not os.path.exists(self._dataPath):
			files = [""]
			pathNotExist = True
		else:
			files = os.listdir(self._dataPath)	
			allDataFiles = filter(self.isDataFiles,files)	
		dataFiles = []
		for f in allDataFiles:	
			s = f.split("_")
			fstep = s[-1]
			if int(fstep) == step:
				dataFiles.append(f)

		nblocks = len(dataFiles)
		blocks = []
		for f in dataFiles:
			s = f.split("_")
			blockIndex = s[1]
			block = self.loadFile(f)
			blocks.append(block)	
		return blocks

	def loadFile(self,fname):
                db = {}
                fn = self._dataPath+'/'+fname
                if not os.path.exists(fn):
                        print "Error: file %s does not exist" % f
                        sys.exit(0)
                f =open(fn)
                header = ["xpos","ypos","u"]
                line = f.readline()
                s= line.split()
		if len(s) > 0:
                	for i,h in enumerate(header):
                        	try:
                                        db.setdefault(h,[]).append(float(s[i]))
           	                except ValueError:
                                        db.setdefault(h,[]).append(float('nan'))
                while line:
                        line = f.readline()
                        s = line.split()
                        if len(s) > 0:
                                for i,h in enumerate(header):
                                        try:
                                                db.setdefault(h,[]).append(float(s[i]))
                                        except ValueError:
                                                db.setdefault(h,[]).append(float('nan'))

                for key in db.keys():
                        db[key] = na.array(db[key])
                return db

if __name__ == '__main__':
	grids = Load(step=98,dataPath="../out").grids
	Plot(grids)
