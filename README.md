# BayesQcd

## rebalance and smear (R&S) code
This code implements a rebalance and smear (R&S) prediction for fake-MET backgrounds by carrying out an event-by-event posterior density maximization. The primary use R&S is to analyze real proton-proton collision data events and transform them into a set that represents the fake-MET background (e.g., QCD, DY->ee/mumu). The method is suitable for dark matter and SUSY searches which require large MET in the analysis search regions.  

## example: make QCD prediction based on simulation
The following example takes as input a provided Delphes QCD (simulated) sample and produces a flat ROOT tree (\textbf{skim}). The output file, when corresponding to real data, can be taken as a replacement for QCD Monte Carlo events. The skim contains both the cloned original events as well as new events which have been rebalanced and smeared - the latter are to be analyzed as if they were MC events in order to create the prediction, and the former, in the case in which they are simulated QCD events, can be compared to the prediction as a cross check (closure test). 

```
input: QCD Delphes event sample ROOT file, which contains an event tree
output: flat skim containing original events and R&S events

==>the \textbf{skim} can be analyzed to produce the prediction and a closure test
```

note: since there are many times more R&S events than original events due to multiple smearing iterations, an alternative use for the simulated R&S events can be as a ready-made background training sample for a multivariate classifier whose goal is to reject fake MET backgrounds. 


## set up code R&S code

Provided that ROOT is installed and pyROOT is enabled, the code can be checked out via

```
git clone https://github.com/sbein/BayesQcd/
cd BayesQcd/
mkdir jobs
mkdir output
mkdir pdfs
```

## run the rebalance and smear code
We can run the code over a QCD file by doing: 

```
python tools/MaximizePosteriorClean.py --quickrun True [--fnamekeyword <path Delphes tree file name, wrap in quotes if using wildcards>]
```

This script takes as input the jet response histograms and prior PDF histograms, which should be remade each time the experiment changes. 
Generate plots overlaying observed and R&S histograms. The option --quickrun species to only analyze the first 5k events. R&S takes 0.1 s/evt, so large data sets and real world problems should be interfaced with a batch system. 

One can open up this file and browse the TTree object (littletree). Note it contains a branch IsRandS, which is a boolean that specifies whether an event is an original un-sullied seed event (false) or a rebalanced and smeared event (true). The latter forms the basis of the prediction, and is what is used in the real data. We can also call this script to draw weighted comparisons between the original and R&S distributions:
```
python tools/DrawAnalyzeClean.py littletree-delphes_qcd_12of80.root
```

This creates both a set of pdfs as well as a root file with canvases. You'll notice your histograms will likely suffer from low statistics, which is why it's good to use the batch system for a large test.

