#!/bin/bash
########################################################################
#
# Shell script to run a suntans test case.
#
########################################################################

SUNTANSHOME=../suntans-particles/main
SUN=$SUNTANSHOME/sun
SUNPLOT=$SUNTANSHOME/sunplot

. $SUNTANSHOME/Makefile.in

maindatadir=rundata
datadir=data

NUMPROCS=$1

if [ -z "$MPIHOME" ] ; then
    EXEC=$SUN
else
    EXEC="$MPIHOME/bin/mpirun -np $NUMPROCS --use-hwthread-cpus $SUN"
fi

if [ ! -d $datadir ] ; then
    cp -r $maindatadir $datadir
    echo Creating grid...
    if [ -z "$TRIANGLEHOME" ] ; then
	echo No triangle libraries installed.  
	echo Copying points.dat, cells.dat, and edges.dat from $maindatadir to $datadir
	$EXEC -g --datadir=$datadir
    else
	$EXEC -t -g --datadir=$datadir
    fi
    echo Placing particles...
    ./make_parts.py
else
    cp $maindatadir/suntans.dat $datadir/.
    cp $maindatadir/lagrange.dat $datadir/.
fi

echo Running suntans...
$EXEC -s -vv --datadir=$datadir

