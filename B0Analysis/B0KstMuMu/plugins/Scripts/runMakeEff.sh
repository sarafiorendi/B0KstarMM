if [ "$#" -ne 1 ] ; then
    echo @@@ Error: file name suffix is missing @@@
else
    source InitAnalysis.sh

    echo @@@ Producing binned efficiency: $1 @@@

    Qsub -l lnxfarm -e -o EffSign_$1.log -N EffSign$1 .././ComputeEfficiency Make $DATADIR/MonteCarlo"$DATAYEAR"/ForEfficiency/B0ToKstMuMu_GEN_NoFilter_MC_NTuples_addGENvars.root $DATADIR/MonteCarlo"$DATAYEAR"/ForEfficiency/B0ToKstMuMu_MC_NTuple_addGENvars.root $DATADIR/MonteCarlo"$DATAYEAR"/SingleCand/singleCand_B0ToKstMuMu_MC_NTuple.root Efficiency4D_B0ToKstMuMu_$1.txt 1

    Qsub -l lnxfarm -e -o EffJPsi_$1.log -N EffJPsi$1 .././ComputeEfficiency Make $DATADIR/MonteCarlo"$DATAYEAR"/ForEfficiency/B0ToJPsiKst_GEN_NoFilter_MC_NTuples_addGENvars.root $DATADIR/MonteCarlo"$DATAYEAR"/ForEfficiency/B0ToJPsiKst_MC_NTuple_addGENvars.root $DATADIR/MonteCarlo"$DATAYEAR"/SingleCand/singleCand_B0ToJPsiKst_MC_NTuple.root Efficiency4D_B0ToJPsiKst_$1.txt 3

    Qsub -l lnxfarm -e -o EffPsi2S_$1.log -N EffPsi2S$1 .././ComputeEfficiency Make $DATADIR/MonteCarlo"$DATAYEAR"/ForEfficiency/B0ToPsi2SKst_GEN_NoFilter_MC_NTuples_addGENvars.root $DATADIR/MonteCarlo"$DATAYEAR"/ForEfficiency/B0ToPsi2SKst_MC_NTuple_addGENvars.root $DATADIR/MonteCarlo"$DATAYEAR"/SingleCand/singleCand_B0ToPsi2SKst_MC_NTuple.root Efficiency4D_B0ToPsi2SKst_$1.txt 5
fi
