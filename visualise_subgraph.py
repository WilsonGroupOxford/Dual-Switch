import sys
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib._color_data as mcd
import matplotlib.colors as colors
import matplotlib.pylab as pylab
import numpy as np

def main():
    inPrefix=sys.argv[1]
    dual, graph, dualColours, graphColours, dualLabels, graphLabels, saveFig, periodic, figSize = visualisationType()
    if(figSize): updateParamsFigSize()
    else: updateParams()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.axis('off')
    if(dual):
        dualCrds, dualCnxs, dualSizes, dualEdges, lattice=readDual(inPrefix,periodic)
        plotDual(dualCrds,dualCnxs,dualSizes,dualEdges,dualColours,dualLabels,fig,ax)
        # setDualAxesLimits(dualCrds,ax,lattice)
    if(graph):
        graphCrds, graphRings, graphSizes, graphImage, lattice=readGraph(inPrefix,periodic)
        # graphCrds[:,[0,1]]=graphCrds[:,[1,0]]
        plotGraph(graphCrds,graphRings,graphSizes,graphImage,graphColours,graphLabels,fig,ax)
        setGraphAxesLimits(graphCrds,graphRings,graphImage,ax,lattice)
    if(saveFig): savePlot(inPrefix)
    displayPlot()
    return

def updateParams():
    params = {'legend.fontsize': 8,
              'font.size' : 8,
              'figure.figsize': (6.5, 6.5),
              'axes.labelsize': 8,
              'axes.titlesize': 8,
              'xtick.labelsize': 8,
              'ytick.labelsize': 8}
    pylab.rcParams.update(params)
    return

def updateParamsFigSize():
    params = {'legend.fontsize': 8,
              'font.size' : 4,
              'figure.figsize': (2.8, 2.8),
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
    graphLabels=False
    save=False
    periodic=False
    figSize=False
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
    if("L" in sys.argv[2]):
        graphLabels=True
    if("s" in sys.argv[2]):
        save=True
    if("p" in sys.argv[2]):
        periodic=True
    if("f" in sys.argv[2]):
        figSize=True
    return dual, graph, dualColours, graphColours, dualLabels, graphLabels, save, periodic, figSize

def readDual(prefix,periodic):
    if(periodic):
        crdFileName=prefix+"_dual_periodic_coordinates.out"
        cnxFileName=prefix+"_dual_periodic_connectivity.out"
        sizeFileName=prefix+"_dual_periodic_size.out"
#        latticeFileName=prefix+"_periodic_lattice_dim.out"
#        lattice=np.genfromtxt(latticeFileName)
        lattice=np.array([1,1])
    else:
        crdFileName=prefix+"_dual_coordinates.out"
        cnxFileName=prefix+"_dual_connectivity.out"
        sizeFileName=prefix+"_dual_size.out"
        lattice=np.array([1,1])
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
    return crds, cnxs, size, edge, lattice

def readGraph(prefix,periodic):
    if(periodic):
        crdFileName=prefix+"_graph_periodic_coordinates.out"
        sizeFileName=prefix+"_graph_periodic_size.out"
        cnxFileName=prefix+"_graph_periodic_rings.out"
        #latticeFileName=prefix+"_periodic_lattice_dim.out"
        #lattice=np.genfromtxt(latticeFileName)
        lattice=np.array([1,1])
    else:
        crdFileName=prefix+"_graph_coordinates.out"
        sizeFileName=prefix+"_graph_size.out"
        cnxFileName=prefix+"_graph_rings.out"
        lattice=np.array([1,1])
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
    else: image=np.ones(len(cnxsList))
    return crds, cnxs, size, image, lattice

def plotDual(crds,cnxs,ringSizes,edges,colourFlag,labelFlag,fig,ax):
    subgraph=7
    for cnxList in cnxs:
        cnx0=cnxList[0]
        for cnx1 in cnxList[1:]:
            if(cnx0<cnx1 and ringSizes[cnx0]==subgraph and ringSizes[cnx1]==subgraph): plt.plot([crds[cnx0][0],crds[cnx1][0]],[crds[cnx0][1],crds[cnx1][1]],color="k",lw=0.5,zorder=2)
    mask=ringSizes==subgraph
    crds=crds[mask,:]
    if(colourFlag):
        ccolours=generateColours(ringSizes,edges,True)
        colours=[]
        for i,c in enumerate(ccolours):
            if mask[i]: colours.append(c)
        plt.scatter(crds[:,0],crds[:,1],c=colours,s=10,zorder=3)
        #plt.scatter(crds[:,0],crds[:,1],edgecolor="k",c=colours,linewidth=1,s=10,zorder=3)
    else: plt.scatter(crds[:,0],crds[:,1],c='k',s=2,zorder=3)
    if(labelFlag):
        for i,crd in enumerate(crds):
            plt.text(crd[0], crd[1],str(i))
    return

def plotGraph(crds,rings,ringSizes,graphImage,colourFlag,labelFlag,fig,ax):
    colours=generateColours(ringSizes)
    polygonCmds=generatePolygonDrawingCommands(4,12);
    colourFilter=generateColourFilter(graphImage)
    nCrds=crds[:,0].size
    # plt.scatter(crds[:,0],crds[:,1],c='b',s=2,zorder=3)
    for i, ring in enumerate(rings):
       ringCrds=np.array([crds[a] for a in ring])
       ringCrds=np.append(ringCrds, [ringCrds[0]], axis=0)
       path=Path(ringCrds, polygonCmds[ringSizes[i]-4])
       patch = patches.PathPatch(path, facecolor=colours[i], lw=0.5, alpha=colourFilter[i])
       ax.add_patch(patch)
       if(labelFlag and graphImage[i]==1):
           ringCom=[np.average(ringCrds[:-1,0])-0.4,np.average(ringCrds[:-1,1])-0.4]
           if(ringCrds[:,0].size-1)==4: ringCom[0]+=0.1
           # print i, ringCrds[:,0].size-1
           # if(i==4):
           #     ringCom[0]+=0.1
           #     # ringCom[1]+=0.15
           # if(i==5):
           #     ringCom[0]-=0.4
           #     ringCom[1]+=0.15
           # if(i==6):
           #     ringCom[0]+=0.3
           #     ringCom[1]-=0.1
           # if(i==9):
           #     # ringCom[0]-=0.1
           #     ringCom[1]+=0.3
           # if(i==10):
           #     ringCom[0]-=0.1
           #     ringCom[1]-=0.1
           # ringCom[0]+=0.8
           # ringCom[1]+=0.5
           plt.text(ringCom[0],ringCom[1],str(ringCrds[:,0].size-1))
       elif(labelFlag):
           print i, ringCrds[:,0].size-1
           ringCom=[np.average(ringCrds[:-1,0])-0.4,np.average(ringCrds[:-1,1])-0.4]
           if(ringCrds[:,0].size-1)==4: ringCom[0]+=0.1
           plt.text(ringCom[0],ringCom[1],str(ringCrds[:,0].size-1),color="dimgrey")

    return

def setDualAxesLimits(crds,ax,lat):
    limLb=np.amin(crds[:,0])-10
    limUb=np.amax(crds[:,0])+10
    sf=lat[1]/lat[0]
    ax.set_xlim(limLb,limUb)
    ax.set_ylim(limLb*sf,limUb*sf)
    return

def setGraphAxesLimits(crds,rings,image,ax,lat):
    n=len(rings)
    ringMask=np.zeros(n,dtype=bool)
    for i, im in enumerate(image):
        if (im==1): ringMask[i]=1
    crdMask=np.zeros(crds[:,0].size,dtype=bool)
    for ring in rings[ringMask]:
        for crd in ring: crdMask[crd]=1
    limLb=np.amin(crds[crdMask,0])-10
    limUb=np.amax(crds[crdMask,0])+10
    sf=lat[1]/lat[0]
    sf=1
    #limLb=-12
    #limUb=24
    ax.set_xlim(limLb,limUb)
    ax.set_ylim(limLb*sf,limUb*sf)
    return

def generateColours(ringSizes,kRings=None,kFlag=False):
    nRings=ringSizes.size
    colours=[]
    colourList=['lightgreen',mcd.CSS4_COLORS['lightblue'],'xkcd:grey',mcd.CSS4_COLORS['salmon'],mcd.CSS4_COLORS['orange'],\
                        mcd.CSS4_COLORS['plum'],mcd.CSS4_COLORS['pink'],"gold","gold"]

    colormapGreens=plt.cm.get_cmap("Greens")
    colormapBlues=plt.cm.get_cmap("Blues")
    colormapGreys=plt.cm.get_cmap("Greys")
    colormapReds=plt.cm.get_cmap("Reds")
    colormapOranges=plt.cm.get_cmap("YlOrBr")
    colormapPurples=plt.cm.get_cmap("PuRd")
    colormapPinks=plt.cm.get_cmap("RdPu")

    green=colormapGreens(100)
    blue=colormapBlues(150)
    grey=colormapGreys(90)
    red=colormapReds(105)
    orange=colormapOranges(100)
    purple=colormapPurples(100)
    pink=colormapPinks(80)

    colourList=[green,blue,grey,red,orange,purple,pink]

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
        if(im!=1): filter[i]=0.2
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
    filename=prefix+".png"
    plt.savefig(filename, dpi=400, bbox_inches="tight")
    return

def displayPlot():
    plt.show()
    return
main()
