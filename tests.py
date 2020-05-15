import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist

colors = ['#DA291CFF','#56A8CBFF','#53A567FF',"#E95C20FF", "#006747FF", "#4F2C1DFF"]

kwargs = dict(alpha=1.0, density=False, stacked=False, hatch='//',linewidth=2, histtype='step')

def getPairs(ref,new):
    appended_data = []
    column_names = ref.get_group(0).columns
    for ev in range(0,len(ref)):
        ref_ev = ref.get_group(ev)
        new_ev = new.get_group(ev)
        distance_array = cdist(ref_ev.iloc[:,:2], new_ev.iloc[:,:2], metric='minkowski', p=2) #returns pairwise distance matrix ( shape=(ref,new), new <= ref)
        # np.set_printoptions(suppress=True,precision=3)
        # print distance_array
        minInRows = np.array(np.argmin(distance_array, axis=1)) #Finds minimum distance between reference object and each new object
        selectedJets = new_ev.iloc[minInRows,:]
        frames = [ref_ev,selectedJets]
        result = pd.concat(frames,axis=1,sort=False)
        appended_data.append(result)
    appended_data = pd.concat(appended_data)
    appended_data.columns = ['ref_' + col_name for col_name in column_names] + ['new_' + col_name for col_name in column_names]
    return appended_data

def getResiduals(ref,new):
    return  (ref - new) / ref

def drawDataFrames(dfs,labels,colors,bins,range,xlabel,ylabel,outpath,outname, log=False):
    plt.clf()
    fig,ax = plt.subplots()
    plt.grid(True)
    for i,df in enumerate(dfs):
        if i == 0:
            ax = df.plot.hist(bins=bins,range=range,edgecolor=colors[i],label=labels[i],**kwargs)
        else:
            df.plot.hist(ax=ax,bins=bins,range=range,edgecolor=colors[i],label=labels[i],**kwargs)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
   
    plt.legend(loc='upper right')
    plt.legend(prop={'size':10}, frameon=False)
    #plt.show()
    axis = plt.gca()
    ymin, ymax = axis.get_ylim()
    ax.set_ylim(0, ymax*2.5)
    ax.autoscale(tight=False)
    if log:
        plt.yscale('log')
        ax.set_ylim(0.1, ymax * 100)
    plt.savefig(outpath+outname+".png")
    plt.savefig(outpath+outname+".pdf")
    plt.clf()

def compareCollections(refCollection,newCollections=[],baseCollectionName='AK4',newCollectionNames=[],outpath="/eos/home-t/thaarres/www/L1T_phase2_PF_jets/"):

    refCollection_event = refCollection.groupby(refCollection.index)
    newCollections_event = []
    matchedJetPairs = []
    nJetsResiduals = []
    labels  = []
    for i,newCollection in enumerate(newCollections):
      newCollection_event = newCollection.groupby(newCollection.index)  
      sigmaN = (refCollection_event.size()- newCollection_event.size())/refCollection_event.size()
      newCollections_event.append(newCollection_event)
      nJetsResiduals.append(sigmaN)
      labels.append('%s <m> = %.3f)'%(newCollectionNames[i],np.mean(sigmaN)))
      
      matchedJetPairs.append(getPairs(refCollection_event,newCollection_event))
    drawDataFrames(nJetsResiduals, labels, colors = colors,bins=100,range=(-.25,.25),xlabel='$N_{jets}^{ref}-N_{jets}^{new}/N_{jets}^{ref}$',ylabel='Number of jets',outpath=outpath, outname="nJets_residuals",log=True)
    drawDataFrames([refCollection_event.size()]+[i.size() for i  in newCollections_event], [baseCollectionName]+newCollectionNames, colors = colors,bins=20,range=(0.,20.),xlabel='$N_{jets}$',ylabel='Number of jets',outpath=outpath, outname="nJets")
    
    columns = refCollection.columns
    fancynames = {'eta':'$\eta$', 'phi':'$\phi$','pt':'$p_T$'}
    for col in columns:
        name = fancynames[col]
        log = False
        x_range=(-3.2,3.2)
        bins = 20
        if col.find('pt')!=-1:
            log = True
            x_range = (0,500)
            bins = 50
        dfs=[]
        labels = []
        for i,matchedJetPair in enumerate(matchedJetPairs):
          residuals=getResiduals(matchedJetPair['ref_'+col],matchedJetPair['new_'+col])
          dfs.append(residuals)
          labels.append('%s (mean = %.3f)'%(newCollectionNames[i],np.mean(residuals)))
        drawDataFrames(dfs,labels,colors = colors,bins=50,range=(-.25,.25),xlabel='%s$^{ref}$-%s$^{new}$/%s$^{ref}$'%(fancynames[col],fancynames[col],fancynames[col]),ylabel='Number of jets',outpath=outpath, outname="%s_residuals"%col,log=True)
        drawDataFrames( [matchedJetPairs[0]['ref_'+col]]+[i['new_'+col] for i  in matchedJetPairs],[baseCollectionName]+newCollectionNames, colors = colors,bins=bins,range=x_range,xlabel=name,ylabel='Number of jets',outpath=outpath, outname="%s"%col,log=log)

    print("Inspect histograms at {}".format(outpath))


