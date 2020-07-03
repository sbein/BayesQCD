from ROOT import *
from utils import *



ftemplates  = TFile('usefulthings/llhd-prior-coarse-1p2width.root')
ftemplates.cd('splines')
ftemplates.ls()

c1 = mkcanvas('c1')

histnames = ['hRTemplate(gPt20.0-25.0, gEta0.0-0.4)', 'hRTemplate(gPt20.0-25.0, gEta2.5-6.0)','hRTemplate(gPt200.0-300.0, gEta0.0-0.4)']

histnames.reverse()
arg = ''
graphs = []
leg = mklegend(x1=.44, y1=.66, x2=.86, y2=.82, color=kWhite)
for ih, hname in enumerate(histnames):
    h = ftemplates.Get('splines/'+hname+'_graph').Clone()
    h.SetTitle('')
    h.GetXaxis().SetTitle('p_{T}^{reco}/p_{T}^{gen}')
    h.GetYaxis().SetTitle('normalized')
    graphs.append(h)
    graphs[-1].SetLineColor(ih+1)
    graphs[-1].SetLineWidth(3)
    graphs[-1].Draw(arg)
    label = hname.split('(')[-1].split(')')[0].replace('gPt', 'p_{T}=').replace(', gEta', ' GeV, |#eta|=')
    leg.AddEntry(graphs[-1], label)
    arg = 'same'

leg.Draw()
c1.Update()
c1.Print('pdfs/PDFs/responses.pdf')
pause()




c1 = mkcanvas('c1')
#c1.SetLogy()
histnames = ['hGenHardMetPtPtB0(ght0.0-10000.0)']
#histnames = ['hGenHardMetDPhiB0(ght0.0-10000.0)']

histnames.reverse()
arg = ''
graphs = []
leg = mklegend(x1=.44, y1=.66, x2=.86, y2=.82, color=kWhite)
for ih, hname in enumerate(histnames):
    h = ftemplates.Get('splines/'+hname+'_graph').Clone()
    g = h.GetHistogram()
    h.SetTitle('')
    h.GetXaxis().SetTitle('generator hard #slash{E}_{T}')
    h.GetYaxis().SetTitle('normalized')
    graphs.append(h)
    graphs[-1].SetLineColor(ih+1)
    graphs[-1].SetLineWidth(3)
    if ih==0: 
        g.GetYaxis().SetRangeUser(0,1.2*g.GetMaximum())
        g.GetXaxis().SetRangeUser(0,100)
        g.Draw()
    arg = 'same'
    graphs[-1].Draw(arg)


leg.Draw()
c1.Update()
c1.Print('pdfs/PDFs/prior.pdf')
#c1.Print('pdfs/PDFs/prior-dphi.pdf')
pause()

