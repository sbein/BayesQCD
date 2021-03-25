# BayesQcd

## purpose of the code
This code implements a rebalance and smear (R&S) prediction for fake-MET backgrounds by carrying out an event-by-event posterior density maximization. The primary use R&S is to analyze real data events and transform them into a set that represents the fake-MET background (e.g., QCD, DY->ee/mumu). This example is based solely on simulated events produced via LO Pythia8 and fed through the Delphes detector simulation program.

The following takes as input a provided Delphes QCD sample and produces a flat ROOT tree (skim). The skims contain both the cloned original events as well as new events which have been rebalanced and smeared. Then in the example, the R&S events are compared with the original events, essentially a closure test of the method. Since there are many times more R&S events than original events due to multiple smearing iterations, the simulated R&S events can also be readily used as a background training sample for a multivariate classifier. 

```
input: QCD Delphes event sample ROOT file, which contains an event tree
output: flat skim containing original events and R&S events
```
## set up code after ensuring an environment in which PYROOT is set

```
git clone https://github.com/sbein/BayesQcd/
cd BayesQcd/
mkdir jobs
mkdir output
mkdir pdfs
```

## run the rebalance and smear code
Straght away, let's run the code over a QCD file by doing: 

```
python tools/MaximizePosteriorClean.py --quickrun True
```

This script takes as input the jet response histograms and prior PDF histograms, which should be remade each time the experiment changes. 
Generate plots overlaying observed and R&S histograms. The option --quickrun species to only analyze the first 5k events. R&S takes .1s/event, so large data sets and real world problems should be interfaced with a batch system. 

One can open up this file and browse the TTree object (littletree). Note it contains a branch IsRandS, which is a boolean that specifies whether an event is an original un-sullied seed event, or a rebalanced and smeared event. The latter forms the basis of the prediction, and is what is used in the real data. We can also call this script to draw weighted comparisons between the original and R&S distributions:
```
python tools/DrawAnalyzeClean.py littletree-delphes_qcd_12of80.root
```

This creates both a set of pdfs as well as a root file with canvases. You'll notice your histograms will likely suffer from low statistics, which is why it's good to use the batch system for a large test.

