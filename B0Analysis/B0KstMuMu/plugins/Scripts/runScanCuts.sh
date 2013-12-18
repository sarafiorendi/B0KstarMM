echo @@@ Scan number: $1 @@@

Qsub -l lnxfarm -e -o Scan$1Cut0.log -N Cut0 .././B0KstMuMuScanCuts 0 3 /nfs/data37/cms/dinardo/Data2012B0KstMuMuResults/MonteCarlo2012/RECOcands/B0ToKstMuMu_MC_NTuple.root /nfs/data37/cms/dinardo/Data2012B0KstMuMuResults/Data2012/B0ToKstMuMu_Data2012ABCD_NTuples.root

Qsub -l lnxfarm -e -o Scan$1Cut1.log -N Cut1 .././B0KstMuMuScanCuts 1 3 /nfs/data37/cms/dinardo/Data2012B0KstMuMuResults/MonteCarlo2012/RECOcands/B0ToKstMuMu_MC_NTuple.root /nfs/data37/cms/dinardo/Data2012B0KstMuMuResults/Data2012/B0ToKstMuMu_Data2012ABCD_NTuples.root

Qsub -l lnxfarm -e -o Scan$1Cut2.log -N Cut2 .././B0KstMuMuScanCuts 2 3 /nfs/data37/cms/dinardo/Data2012B0KstMuMuResults/MonteCarlo2012/RECOcands/B0ToKstMuMu_MC_NTuple.root /nfs/data37/cms/dinardo/Data2012B0KstMuMuResults/Data2012/B0ToKstMuMu_Data2012ABCD_NTuples.root

Qsub -l lnxfarm -e -o Scan$1Cut3.log -N Cut3 .././B0KstMuMuScanCuts 3 3 /nfs/data37/cms/dinardo/Data2012B0KstMuMuResults/MonteCarlo2012/RECOcands/B0ToKstMuMu_MC_NTuple.root /nfs/data37/cms/dinardo/Data2012B0KstMuMuResults/Data2012/B0ToKstMuMu_Data2012ABCD_NTuples.root

Qsub -l lnxfarm -e -o Scan$1Cut4.log -N Cut4 .././B0KstMuMuScanCuts 4 3 /nfs/data37/cms/dinardo/Data2012B0KstMuMuResults/MonteCarlo2012/RECOcands/B0ToKstMuMu_MC_NTuple.root /nfs/data37/cms/dinardo/Data2012B0KstMuMuResults/Data2012/B0ToKstMuMu_Data2012ABCD_NTuples.root

Qsub -l lnxfarm -e -o Scan$1Cut5.log -N Cut5 .././B0KstMuMuScanCuts 5 3 /nfs/data37/cms/dinardo/Data2012B0KstMuMuResults/MonteCarlo2012/RECOcands/B0ToKstMuMu_MC_NTuple.root /nfs/data37/cms/dinardo/Data2012B0KstMuMuResults/Data2012/B0ToKstMuMu_Data2012ABCD_NTuples.root
