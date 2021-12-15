import os, sys
from ROOT import *
from utils import *
from glob import glob
gROOT.SetBatch(1)
gStyle.SetOptStat(0)

'''
e.g.,
python tools/DrawAnalyze.py "output/smallchunks/littletree*GJet*.root"

python tools/DrawAnalyze.py "output/smallchunks/tree*qcd*.root"

#follow up with:
#hadd weightedHists_littletreeAll_YesAcme.root weightedHists_littletreeqcd_YesAcme.root weightedHists_littletreeWZ_YesAcme.root
'''

try: fileskey = sys.argv[1]
except:
    print 'please enter path to input files, e.g.'
    print 'python tools/DrawAnalyze.py "output/smallchunks/littletree*GJet*.root"'

print 'fileskey', fileskey


universalconstraint = '1 == 1'#' abs(HardMetMinusMet)<80 && mva_Photons1Et>80'
acmestr = 'NoAcme' if 'NoAcme' in fileskey else 'YesAcme'

fins = glob(fileskey)
ccounter = TChain('tcounter')
for fname in fins: ccounter.Add(fname)
nev_total = ccounter.GetEntries()


chain = TChain('littletree')
print 'fileskey', fileskey
for fname in fins: chain.Add(fname)
chain.Show(0)
print 'nevents =', chain.GetEntries(), nev_total

evtweight = 'CrossSectionPb/'+str(nev_total)


plotBundle = {}

plotBundle['Baseline_NJets'] = ['NJets>>hadc(9,1,10)','HardMetPt>120 && BTags==0',True]
plotBundle['Baseline_HardMet'] = ['HardMetPt>>hadc(25,0,500)','1==1 && BTags==0',True]
plotBundle['Baseline_HT'] = ['HT>>hadc(28,100,2900)','HardMetPt>120 && BTags==0',True]
plotBundle['Baseline_MinDPhiJetsHardMet'] = ['MinDPhiJetsHardMet>>hadc(28,0.0,3.2)','HardMetPt>120 && BTags==0',True]
    

print 'fileskey', fileskey
infilekey = fileskey.split('/')[-1].replace('*','').replace('.root','')

fnew = TFile('weightedHists_'+infilekey+'_'+acmestr+'.root', 'recreate')
print 'will make file', fnew.GetName()
c1 = mkcanvas()
c2 = mkcanvas('c2')

for key in plotBundle:
    drawarg, constraint, usernsvalue = plotBundle[key]
    obsweight = evtweight+'*('+constraint + ' && '+ universalconstraint + ' && IsRandS==0)'
    #puWeight
    print 'drawing', drawarg, ', with constraint:', obsweight
    chain.Draw(drawarg,obsweight, 'e')
    hobs = chain.GetHistogram().Clone(key+'_obs')
    hobs.GetYaxis().SetRangeUser(0.0001,200*hobs.GetMaximum())

    c1.cd()

    drawarg = drawarg
    randsconstraint = constraint
    methweight = evtweight+'/NSmearsPerEvent*('+ randsconstraint + ' && '+universalconstraint+ ' && IsRandS==1)'
    #puWeight
    print 'drawing', drawarg, ', with constraint:', methweight
    chain.Draw(drawarg, methweight, 'e')
    hrands = chain.GetHistogram().Clone(key+'_rands') 
    hrands.GetYaxis().SetRangeUser(0.1,2000*hrands.GetMaximum())
    if 'ZGG' in fileskey: histoStyler(hrands, kViolet+1)
    else: histoStyler(hrands, kAzure-8)
    hrands.SetFillColor(hrands.GetLineColor())
    hrands.SetFillStyle(1001)

    leg = mklegend(x1=.45, y1=.57, x2=.95, y2=.74, color=kWhite)
    if 'ZGG' in fileskey: hobs.SetTitle('Summer16 ZGG')
    else: hobs.SetTitle('observed')
    hobs.GetXaxis().SetTitle(key.split('_')[-1])
    hrands.GetXaxis().SetTitle(key.split('_')[-1])
    hrands.SetTitle('rebalance and smeared')
    hratio, hmethodsyst = FabDrawSystyRatio(c1,leg,hobs,[hrands],datamc='MC',lumi='n/a', title = '', LinearScale=False, fractionthing='data / prediction')
    c1.Update()
    c1.Write('c_'+key)
    hrands.Write('h'+hrands.GetName())
    hobs.Write('h'+hobs.GetName())
    c1.Print('pdfs/ClosureSkims/'+key+'.pdf')        


print 'just created', fnew.GetName()
fnew.Close()



