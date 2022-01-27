from ROOT import *
from utils import *
gROOT.SetBatch(1)
gStyle.SetOptStat(0)

#ftemplates  = TFile('usefulthings/llhd-prior-coarse-1p2width.root')
ftemplates = TFile('usefulthings/llhd-prior-coarse-realisticHcal.root')
ftemplates.cd('splines')
ftemplates.ls()

c1 = mkcanvas('c1')

histnames = ['hRTemplate(gPt20.0-25.0, gEta0.0-0.4)', 'hRTemplate(gPt20.0-25.0, gEta2.5-6.0)','hRTemplate(gPt200.0-300.0, gEta0.0-0.4)']

histnames.reverse()
arg = 'hist l'
graphs = []
leg = mklegend(x1=.41, y1=.62, x2=.88, y2=.85, color=kWhite)
for ih, hname in enumerate(histnames):
    h = ftemplates.Get(hname).Clone()
    h.Rebin(2)
    h.Scale(1.0/h.Integral())
    h.SetTitle('')
    h.GetXaxis().SetTitle('p_{T}^{reco}/p_{T}^{gen}')
    h.GetYaxis().SetTitle('unit normalized')
    h.Smooth(100)
    histoStyler(h, ih+1)
    h.GetYaxis().SetLabelSize(0.05)
    h.GetYaxis().SetTitleOffset(1.15)    
    graphs.append(h)
    graphs[-1].SetLineColor(ih+1)
    graphs[-1].SetLineWidth(3)
    graphs[-1].Draw(arg)
    label = hname.split('(')[-1].split(')')[0].replace('gPt', 'p_{T}=').replace(', gEta', ' GeV, |#eta|=')
    label = label.replace('.0','')
    leg.AddEntry(graphs[-1], label)
    arg = 'same hist l'

leg.Draw()
c1.Update()
c1.Print('pdfs/PDFs/responses.pdf')
#pause()





#c1.SetLogy()
histnames = ['hGenHardMetPtPtB0(ght0.0-10000.0)', 'hGenHardMetPtPtB1(ght0.0-10000.0)','hGenHardMetPtPtB2(ght0.0-10000.0)']
histnames = ['hGenHardMetDPhiB0(ght0.0-10000.0)','hGenHardMetDPhiB1(ght0.0-10000.0)','hGenHardMetDPhiB2(ght0.0-10000.0)']

#histnames.reverse()
arg = 'lhist'
graphs = []

c1 = mkcanvas('c1')
leg = mklegend(x1=.34, y1=.64, x2=.86, y2=.83, color=kWhite)
for ih, hname in enumerate(histnames):
    #h = ftemplates.Get('splines/'+hname+'_graph').Clone()
    print ih, 'doin', hname
    h = ftemplates.Get(hname).Clone()
    h.Scale(1.0/h.Integral())
    #tSmoother = TGraphSmooth()
    #for i in range(100): tSmoother.Smoothin(g)
    #h = g.GetHistogram()
    histoStyler(h, ih+1)
    #histoStyler(g, 1)
    h.SetTitle('')
    
    h.GetYaxis().SetTitle('unit normalized')
    h.GetYaxis().SetLabelSize(0.045)
    h.GetYaxis().SetTitleOffset(1.12)
    leglabel = 'n(b-jets)='+str(ih)
    if ih==len(histnames)-1:
        leglabel = leglabel.replace('=','#geq')
        histoStyler(h, ih+2)
    #h.GetYaxis().SetRangeUser(0,1.2*h.GetMaximum())
    
    if 'MetPt' in hname:# ih==0:
        h.Rebin(2) 
        h.GetXaxis().SetRangeUser(0,100)
        h.GetXaxis().SetTitle('#slash{p}_{T}^{hard} [GeV]')
    else:
        #h.Smooth(100)
        h.Rebin(3)
        h.GetXaxis().SetTitle('#Delta#phi(#slash{#vec{p}}_{T}^{hard}, #vec{p}_{T}^{j1})')
        a = 2
    h.GetYaxis().SetRangeUser(0,1.2*h.GetMaximum())        
    h.Smooth(100)
    graphs.append(h)
    #graphs[-1].SetLineColor(1)
    graphs[-1].SetLineWidth(3)
    leg.AddEntry(graphs[-1], leglabel)
    graphs[-1].Draw(arg)
    arg = 'lhist same'    
    #g.Draw('cl same')


leg.Draw()
c1.Update()
if 'DPhi' in histnames[0]: c1.Print('pdfs/PDFs/prior_dphi'+'.pdf')
else: c1.Print('pdfs/PDFs/prior_metpt'+'.pdf')
#c1.Print('pdfs/PDFs/prior-dphi.pdf')
#pause()
