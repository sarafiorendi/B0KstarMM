if [ "$#" -ne 1 ] ; then
    echo @@@ Error: Analysis path is missing @@@
else
    export ANALYPATH=$1

    export DATAYEAR=2012
    export DATADIR=/nfs/data37/cms/dinardo/Data"$DATAYEAR"B0KstMuMuResults

    unset DISPLAY

    echo @@@ Analysis environment variable: $ANALYPATH @@@
    echo @@@ Directory with data: $DATADIR @@@
    echo @@@ Data year: $DATAYEAR @@@
    echo @@@ Unset DISPLAY @@@
fi
