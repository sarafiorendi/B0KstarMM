#!/usr/bin/perl


### GEN MC ###
#$dirName = "/nfs/data36/cms/dinardo/B0ToKstMuMu_GEN_NoFilter_MC_NTuples/" ;
#$dirName = "/nfs/data36/cms/dinardo/B0ToJPsiKst_GEN_NoFilter_MC_NTuples/" ;
#$dirName = "/nfs/data36/cms/dinardo/B0ToPsi2SKst_GEN_NoFilter_MC_NTuples/" ;

### RECO MC ###
#$dirName = "/bestman/storage/cms/store/user/dinardo/GridB0KstMuMu_B0ToKstMuMu/" ;
#$dirName = "/bestman/storage/cms/store/user/dinardo/GridB0KstMuMu_B0ToJPsiKst/" ;
#$dirName = "/bestman/storage/cms/store/user/dinardo/GridB0KstMuMu_B0ToPsi2SKst/" ;
#$dirName = "/bestman/storage/cms/store/user/dinardo/GridB0KstMuMu_BsToKstMuMu/" ;

#$dirName = "/bestman/storage/cms/store/user/dinardo/GridB0KstMuMu_B0ToPsiMuMu/" ;
#$dirName = "/bestman/storage/cms/store/user/dinardo/GridB0KstMuMu_BpToPsiMuMu/" ;
#$dirName = "/bestman/storage/cms/store/user/dinardo/GridB0KstMuMu_BsToPsiMuMu/" ;
#$dirName = "/bestman/storage/cms/store/user/dinardo/GridB0KstMuMu_LambdaBToPsiMuMu/" ;

### Data ###
$dirName = "/bestman/storage/cms/store/user/dinardo/GridB0KstMuMu_Data2012A/" ;


$listcmd = "ls " . $dirName ;
@list = `$listcmd` ;

$nFiles = @list ;
$split = 700 ;
$nOut = int($nFiles / $split) + 1 ;

print "Total number of files: " . $nFiles . "\n" ;
print "Split into: " . $nOut . " files\n" ;

for ($count = 0; $count < $nOut; $count++)
{


### GEN MC ###
#    $fileOut = "B0ToKstMuMu_GEN_NoFilter_MC_NTuples" ;
#    $fileOut = "B0ToJPsiKst_GEN_NoFilter_MC_NTuples" ;
#    $fileOut = "B0ToPsi2SKst_GEN_NoFilter_MC_NTuples" ;

### RECO MC ###
#    $fileOut = "B0ToKstMuMu_MC_NTuple" ;
#    $fileOut = "B0ToJPsiKst_MC_NTuple" ;
#    $fileOut = "B0ToPsi2SKst_MC_NTuple" ;
#    $fileOut = "BsToKstMuMu_MC_NTuple" ;

#    $fileOut = "B0ToPsiMuMu_MC_NTuple" ;
#    $fileOut = "BpToPsiMuMu_MC_NTuple" ;
#    $fileOut = "BsToPsiMuMu_MC_NTuple" ;
#    $fileOut = "LambdaBToPsiMuMu_MC_NTuple" ;

### Data ###
    $fileOut = "B0ToKstMuMu_Data2012A_NTuples" ;


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
