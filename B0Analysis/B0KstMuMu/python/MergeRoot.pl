#!/usr/bin/perl


### RECO Central MC ###
#$dirName = "/nfs/data36/cms/dinardo/B0ToPsiMuMu_MC_NTuples/" ;
#$dirName = "/nfs/data36/cms/dinardo/BpToPsiMuMu_MC_NTuples/" ;
#$dirName = "/nfs/data36/cms/dinardo/BsToPsiMuMu_MC_NTuples/" ;
#$dirName = "/nfs/data36/cms/dinardo/LambdaBToPsiMuMu_MC_NTuples/" ;
### GEN MC ###
#$dirName = "/nfs/data36/cms/dinardo/B0ToKstMuMu_GEN_Filter_MC_NTuples/" ;
#$dirName = "/nfs/data36/cms/dinardo/B0ToKstMuMu_GEN_NoFilter_01_MC_NTuples/" ;
#$dirName = "/nfs/data36/cms/dinardo/B0ToJPsiKst_GEN_Filter_MC_NTuples/" ;
#$dirName = "/nfs/data36/cms/dinardo/B0ToJPsiKst_GEN_NoFilter_MC_NTuples/" ;
#$dirName = "/nfs/data36/cms/dinardo/B0ToPsi2SKst_GEN_Filter_MC_NTuples/" ;
#$dirName = "/nfs/data36/cms/dinardo/B0ToPsi2SKst_GEN_NoFilter_MC_NTuples/" ;
### Data ###
$dirName = "/bestman/storage/cms/store/user/dinardo/GridB0KstMuMu_PromptReco_v4A_67/" ;


$listcmd = "ls " . $dirName ;
@list = `$listcmd` ;

$nFiles = @list ;
$split = 700 ;
$nOut = int($nFiles / $split) + 1 ;

print "Total number of files: " . $nFiles . "\n" ;
print "Split into: " . $nOut . " files\n" ;

for ($count = 0; $count < $nOut; $count++)
{


### RECO Central MC ###
#    $fileOut = "B0ToPsiMuMu_MC_NTuples_Merged" ;
#    $fileOut = "BpToPsiMuMu_MC_NTuples_Merged" ;
#    $fileOut = "BsToPsiMuMu_MC_NTuples_Merged" ;
#    $fileOut = "LambdaBToPsiMuMu_MC_NTuples_Merged" ;
### RECO MC ###
#    $fileOut = "B0ToKstMuMu_BrMC_NTuples_Merged" ;
#    $fileOut = "B0ToPsi2SKst_BrMC_NTuples_Merged" ;
#    $fileOut = "B0ToJPsiKst_BrMC_NTuples_Merged" ;
### GEN MC ###
#    $fileOut = "B0ToKstMuMu_GEN_Filter_MC_NTuples" ;
#    $fileOut = "B0ToKstMuMu_GEN_NoFilter_01_MC_NTuples" ;
#    $fileOut = "B0ToJPsiKst_GEN_Filter_MC_NTuples" ;
#    $fileOut = "B0ToJPsiKst_GEN_NoFilter_MC_NTuples" ;
#    $fileOut = "B0ToPsi2SKst_GEN_Filter_MC_NTuples" ;
#    $fileOut = "B0ToPsi2SKst_GEN_NoFilter_MC_NTuples" ;
### Data ###
    $fileOut = "B0ToKstMuMu_DataPRv4_NTuples_Merged" ;


    $cmd = "hadd " . $fileOut . "_" . $count . ".root " ;
    $Start = $split * $count ;
    $End = $Start + $split - 1 ;
    if ($count == ($nOut-1))
    {
	$End = $Start + ($nFiles % $split) - 1 ;
    }

    print "\nStarting file number: " . $Start . "\n" ;
    print "Ending file number: " . $End . "\n" ;

    foreach $file (@list[$Start..$End])
    {
	chomp $file ;
	$cmd .= $dirName . $file . " " ;
    }

    $fileOut .= "_" . $count . ".sh" ;

    print "\nFile output name: " . $fileOut . "\n" ;

    open(OUT, ">" . $fileOut) ;
    print OUT "$cmd" ;
    print OUT "\n" ;
    close(OUT) ;
}
