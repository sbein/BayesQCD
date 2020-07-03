from ROOT import *
from utils import *
from array import array

files2compare = ['usefulthings/ResponseFunctionsMC17AllFilters_deepCsv.root',\
                  'usefulthings/llhd-prior-coarse-1p2width.root',\
                  'usefulthings/llhd-prior-coarse-widenedecalhcal_step4.root',\
                  'usefulthings/llhd-prior-coarse-widenedecalhcal_step6.root']

colors = [kRed, kBlue, kGreen, kCyan]
fgather = TFile(files2compare[0])
fgather.ls()
etaranges = []
keys = fgather.GetListOfKeys()
for key in keys:
    name = key.GetName()
    if not 'hRTemplate' in name: continue
    if ')B' in name: continue
    etarange = name.split()[1]
    if not etarange in etaranges: etaranges.append(etarange)
    
fgather.Close()

print 'eta ranges:', etaranges

#c1.SetLogy()
for etarange in etaranges:
    graphs = []
    c1 = mkcanvas('c1')
    for fname in files2compare:###probably need to switch order between file/eta and so forth
        f = TFile(fname)
        names = []
        keys = f.GetListOfKeys()
        for key in keys:
            name = key.GetName()
            if not 'hRTemplate' in name: continue
            if ')B' in name: continue
            names.append(name)        
        
        x, y = array( 'd' ), array( 'd' )        
        for name in names:
            if not etarange in name: continue
            ptrange = name.split(',')[0].split('gPt')[-1].split('-')
            ptrange = 0.5*(2*float(ptrange[0])+0*float(ptrange[1]))
            x.append(ptrange)
            hist = f.Get(name)
            y.append(hist.GetRMS())
        gr = TGraph( len(x)-1, x, y )
        gr.SetTitle(etarange)
        #gr.SetDirectory(0)
        graphs.append(gr)
        f.Close()
    tl = mklegend(x1=.22, y1=.7, x2=.79, y2=.82, color=kWhite)
    arg = ''
    for ig, graph in enumerate(graphs):
        graph.SetLineColor(colors[ig])
        graph.Draw(arg)
        graph.SetLineWidth(2)
        tl.AddEntry(graph,files2compare[ig].split('/')[-1], 'l')
        arg = 'same'
    tl.Draw()
    c1.Update()
    pause()

