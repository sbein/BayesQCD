python tools/SubmitJobs_condor.py --analyzer tools/MaximizePosterior.py --fnamekeyword "/nfs/dust/cms/user/beinsam/RebalanceAndSmear/CMSSW_10_1_0/src/BayesQcd/../SampleProduction/delphes/output_ak5/delphes_gjet*.root" --quickrun True &

python tools/SubmitJobs_condor.py --analyzer tools/MaximizePosterior.py --fnamekeyword "/nfs/dust/cms/user/beinsam/RebalanceAndSmear/CMSSW_10_1_0/src/BayesQcd/../SampleProduction/delphes/delphes_qcdextra*.root" --quickrun True &

python tools/SubmitJobs_condor.py --analyzer tools/MaximizePosterior.py --fnamekeyword "/nfs/dust/cms/user/beinsam/RebalanceAndSmear/CMSSW_10_1_0/src/BayesQcd/../SampleProduction/delphes/delphes_T2qq.root" &

#python tools/SubmitJobs_condor.py --analyzer tools/MaximizePosterior.py --fnamekeyword "/nfs/dust/cms/user/beinsam/RebalanceAndSmear/CMSSW_10_1_0/src/BayesQcd/../SampleProduction/delphes/delphes_T2qqGMSB*.root" &


##large number of qcd:
python tools/SubmitJobs_condor.py --analyzer tools/MaximizePosterior.py --fnamekeyword "/nfs/dust/cms/user/beinsam/RebalanceAndSmear/CMSSW_10_1_0/src/SampleProduction/delphes/rootfiles/*.root" --quickrun True &



