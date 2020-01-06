import pandas as pd
import uproot
import numpy as np

def loadFromFile(filename,writeToh5=False):
    f = uproot.open(filename)
    puppi_tree = 'Events/l1tPFCandidates_l1pfCandidates_Puppi_RESP./l1tPFCandidates_l1pfCandidates_Puppi_RESP.obj/'
    p4_tree = puppi_tree + 'l1tPFCandidates_l1pfCandidates_Puppi_RESP.obj.m_state.p4Polar_.fCoordinates'
    pt = f[p4_tree + '.fPt'].array()
    eta = f[p4_tree + '.fEta'].array()
    phi = f[p4_tree + '.fPhi'].array()

    index = np.array([i for i, ev in enumerate(pt) for j in range(len(ev))]).flatten()
    events = pd.DataFrame({'i' : index, 'pt' : pt.flatten(), 'eta' : eta.flatten(), 'phi' : phi.flatten()})
    events.set_index('i', inplace=True, drop=True)

    ak4_tree = 'Events/recoPFJets_ak4PFL1Puppi__RESP./recoPFJets_ak4PFL1Puppi__RESP.obj'
    p4_tree = ak4_tree + '/recoPFJets_ak4PFL1Puppi__RESP.obj.m_state.p4Polar_.fCoordinates'
    pt = f[p4_tree + '.fPt'].array()
    eta = f[p4_tree + '.fEta'].array()
    phi = f[p4_tree + '.fPhi'].array()
    index = np.array([i for i, ev in enumerate(pt) for j in range(len(ev))]).flatten()
    jets = pd.DataFrame({'i' : index, 'pt' : pt.flatten(), 'eta' : eta.flatten(),      'phi' : phi.flatten()})
    jets.set_index('i', inplace=True, drop=True)
    if writeToh5:
        with pd.HDFStore(filename.replace('.root','.h5'),mode='w') as f:
            f['events'] = events
            f['jets']   = jets
    return events, jets

def plotJetCompare(particles, jets1, jets2, label1='Jets 1', label2='Jets 2', radius=0.4):
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.patches import Circle
    scale = 1.
    plt.figure()
    ph = plt.scatter(particles['eta'], particles['phi'], s = particles['pt'] / scale, label='Particles')
    ax = plt.gca()
    for iJ, j in enumerate(jets1.iterrows()):
         j = j[1]
         plt.plot(j['eta'], j['phi'], marker='x', color='blue')
         circle = Circle((j['eta'], j['phi']), radius=radius, edgecolor='blue', fill=False)
         ax.add_patch(circle)
    for iJ, j in enumerate(jets2.iterrows()):
         j = j[1]
         plt.plot(j['eta'], j['phi'], marker='x', color='orange')
         circle = Circle((j['eta'], j['phi']), radius=radius, edgecolor='orange', fill=False, linestyle='--')
         ax.add_patch(circle)
    patch1 = mpatches.Patch(color='blue', label=label1)
    patch2 = mpatches.Patch(color='orange', label=label2)
    plt.legend(handles=[ph, patch1, patch2])
    plt.axis('equal')
    plt.xlabel('eta')
    plt.ylabel('phi')
