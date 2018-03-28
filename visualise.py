import sys
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib._color_data as mcd
import matplotlib.pylab as pylab
import numpy as np

def main():
    inPrefix=sys.argv[1]
    dual, graph, dualColours, graphColours, dualLabels, saveFig, periodic = visualisationType()
    updateParams()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.axis('off')
    if(dual):
        dualCrds, dualCnxs, dualSizes, dualEdges=readDual(inPrefix,periodic)
        plotDual(dualCrds,dualCnxs,dualSizes,dualEdges,dualColours,dualLabels,fig,ax)
        setDualAxesLimits(dualCrds,ax)
    if(graph):
        graphCrds, graphRings, graphSizes, graphImage=readGraph(inPrefix,periodic)
        plotGraph(graphCrds,graphRings,graphSizes,graphImage,graphColours,fig,ax)
        #setGraphAxesLimits(graphCrds,graphRings,graphImage,ax)
    if(saveFig): savePlot(inPrefix)
    displayPlot()
    return

def updateParams():
    params = {'legend.fontsize': 8,
              'figure.figsize': (6.5, 6.5),
              'axes.labelsize': 8,
              'axes.titlesize': 8,
              'xtick.labelsize': 8,
              'ytick.labelsize': 8}
    pylab.rcParams.update(params)
    return

def visualisationType():
    dual=False
    graph=False
    dualColours=False
    graphColours=False
    dualLabels=False
    save=False
    periodic=False
    if("d" in sys.argv[2]):
        dual=True
    if("g" in sys.argv[2]):
        graph=True
    if("c" in sys.argv[2]):
        dualColours=True
    if("C" in sys.argv[2]):
        graphColours=True
    if("l" in sys.argv[2]):
        dualLabels=True
    if("s" in sys.argv[2]):
        save=True
    if("p" in sys.argv[2]):
        periodic=True
    return dual, graph, dualColours, graphColours, dualLabels, save, periodic

def readDual(prefix,periodic):
    if(periodic):
        crdFileName=prefix+"_dual_periodic_coordinates.out"
        cnxFileName=prefix+"_dual_periodic_connectivity.out"
    	sizeFileName=prefix+"_dual_periodic_size.out"
    else:
        crdFileName=prefix+"_dual_coordinates.out"
        cnxFileName=prefix+"_dual_connectivity.out"
    	sizeFileName=prefix+"_dual_size.out"
    crds=np.genfromtxt(crdFileName)
    sizeAndEdge=np.genfromtxt(sizeFileName,dtype=int)
    size=sizeAndEdge[:,0]
    edge=sizeAndEdge[:,1]
    cnxFile = open(cnxFileName, "r")
    cnxsList = []
    for line in cnxFile:
        cnx = [int(x) for x in line.split()]
        cnxsList.append(cnx)
    cnxFile.close()
    cnxs=np.array(cnxsList)
    return crds, cnxs, size, edge

def readGraph(prefix,periodic):
    if(periodic):
        crdFileName=prefix+"_graph_periodic_coordinates.out"
        sizeFileName=prefix+"_graph_periodic_size.out"
        cnxFileName=prefix+"_graph_periodic_rings.out"
    else:
        crdFileName=prefix+"_graph_coordinates.out"
        sizeFileName=prefix+"_graph_size.out"
        cnxFileName=prefix+"_graph_rings.out"
    crds=np.genfromtxt(crdFileName)
    size=np.genfromtxt(sizeFileName,usecols=0,dtype=int)
    cnxFile = open(cnxFileName, "r")
    cnxsList = []
    for line in cnxFile:
        cnx = [int(x) for x in line.split()]
        cnxsList.append(cnx)
    cnxFile.close()
    cnxs=np.array(cnxsList)
    if(periodic):
        size=np.tile(size,9)
        image=np.genfromtxt(sizeFileName,usecols=1,dtype=int)
    else: image=np.array([4]*len(cnxsList))
    return crds, cnxs, size, image

def plotDual(crds,cnxs,ringSizes,edges,colourFlag,labelFlag,fig,ax):
    for cnxList in cnxs:
        cnx0=cnxList[0]
        for cnx1 in cnxList[1:]:
            if(cnx0<cnx1): plt.plot([crds[cnx0][0],crds[cnx1][0]],[crds[cnx0][1],crds[cnx1][1]],color="k",lw=0.5,zorder=2)
    if(colourFlag):
    	colours=generateColours(ringSizes,edges,True)
	plt.scatter(crds[:,0],crds[:,1],c=colours,s=2,zorder=3)
    else: plt.scatter(crds[:,0],crds[:,1],c='k',s=2,zorder=3)
    if(labelFlag):
        for i,crd in enumerate(crds):
            plt.text(crd[0], crd[1],str(i))
    return

def plotGraph(crds,rings,ringSizes,graphImage,colourFlag,fig,ax):
    #colours=generateColours(ringSizes)
    #polygonCmds=generatePolygonDrawingCommands(4,12);
    #colourFilter=generateColourFilter(graphImage)
    #nCrds=crds[:,0].size
    plt.scatter(crds[:,0],crds[:,1],c='b',s=2,zorder=3)
    #for i, ring in enumerate(rings):
    #   ringCrds=np.array([crds[a] for a in ring])
    #   ringCrds=np.append(ringCrds, [ringCrds[0]], axis=0)
    #   path=Path(ringCrds, polygonCmds[ringSizes[i]-4])
    #   patch = patches.PathPatch(path, facecolor=colours[i], lw=1, alpha=colourFilter[i])
    #   ax.add_patch(patch)
    return

def setDualAxesLimits(crds,ax):
    limLb=np.amin(crds[:,0])-10
    limUb=np.amax(crds[:,0])+10
    ax.set_xlim(limLb,limUb)
    ax.set_ylim(limLb,limUb)
    return

def setGraphAxesLimits(crds,rings,image,ax):
    n=len(rings)
    ringMask=np.zeros(n,dtype=bool)
    for i, im in enumerate(image):
        if (im==4): ringMask[i]=1
    crdMask=np.zeros(crds[:,0].size,dtype=bool)
    for ring in rings[ringMask]:
        for crd in ring: crdMask[crd]=1
    limLb=np.amin(crds[crdMask,0])-25
    limUb=np.amax(crds[crdMask,0])+25
    ax.set_xlim(limLb,limUb)
    ax.set_ylim(limLb,limUb)
    return;

def generateColours(ringSizes,kRings=None,kFlag=False):
    nRings=ringSizes.size
    colours=[]
    colourList=['lightgreen',mcd.CSS4_COLORS['lightblue'],'xkcd:grey',mcd.CSS4_COLORS['salmon'],mcd.CSS4_COLORS['orange'],\
                        mcd.CSS4_COLORS['plum'],mcd.CSS4_COLORS['pink'],"gold","gold"]
    if(kFlag):
        for i in range(nRings):
	    if(int(ringSizes[i])<4 or kRings[i]==1): colours.append('black')
	    else: colours.append(colourList[int(ringSizes[i])-4])
    else:
        for i in range(nRings):
	    if(int(ringSizes[i])<4): colours.append('black')
	    else: colours.append(colourList[int(ringSizes[i])-4])
    return colours

def generateColourFilter(image):
    n=image.size
    filter=np.ones(n,dtype=float)
    for i, im in enumerate(image):
        if(im!=4): filter[i]=0.5
    return filter

def generatePolygonDrawingCommands(min, max):
    polygonCode = [1, 79]
    for i in range(1, min):
        polygonCode.insert(1, 2)
    allPolyCodes = []
    allPolyCodes.append(polygonCode[:])
    for i in range(min, max):
        polygonCode.insert(1, 2)
        allPolyCodes.append(polygonCode[:])
    return allPolyCodes

def savePlot(prefix):
    filename=prefix+".pdf"
    plt.savefig(filename, dpi=400, bbox_inches="tight")
    return

def displayPlot():
    plt.show()
    return
main()
