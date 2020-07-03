from ROOT import *
from utils import *
gStyle.SetOptStat(0)


batchmode = True

if batchmode: gROOT.SetBatch(1)

fname1 = 'ma5plotsT5Pythia1800.root'
fname2 = 'ma5plotsT5Mad1800.root'


f1 = TFile(fname1)
f2 = TFile(fname2)


f1.ls()

histnames = ['hNLepton','hPtGenLeptons','hEtaGenLeptons','hPtRecoLeptons','hPtRecoIsoLeptons']


c1 = mkcanvas('c1')
for histname in histnames:
    h1 = f1.Get(histname)
    h2 = f2.Get(histname) 
    histoStyler(h1,kBlack)
    histoStyler(h2,kTeal-1)
    h1.SetTitle('pythia')
    h2.SetTitle('mg')    
    h1.Scale(1./h1.Integral(-1,999))
    h2.Scale(1./h2.Integral(-1,999))    
    leg = mklegend(x1=.62, y1=.66, x2=.89, y2=.82, color=kWhite)
    hratio = FabDraw(c1,leg,h1,[h2],datamc='MC',lumi='none', title = '', LinearScale=False, fractionthing= 'pyth/mad')
    hratio.GetYaxis().SetRangeUser(0,2)
    hratio.GetXaxis().SetTitle(histname)    
    c1.Update()
    c1.Print('pdfs/T5tests/'+histname+'.pdf')
    if not batchmode: pause()




