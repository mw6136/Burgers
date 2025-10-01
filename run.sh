#!/bin/bash

RUNNAME="example_name"
TIMEINT=1    # 0 for explicit, 1 for implicit
SPATINT=1     # 0 for centered, 1 for upwind
DIFFCOEF=.01  # Diffusion coefficient
NX1=128        # Number of grid points

MP4SAVE=1     # 0 does not convert an mp4 movie of simulation, 1 does convert
MOVIEFPS=20


# ========================================================================================================= #


if [ ! -d $RUNNAME ]
then
    mkdir $RUNNAME
fi

cd $RUNNAME

if [ ! -d outputs ]
then
    mkdir outputs
fi

clang++ -I/opt/homebrew/include -L/opt/homebrew/lib -larmadillo -std=c++11 ../src/main.cpp -o burger

if [ -f "burger" ]
then
    ./burger ${TIMEINT} ${SPATINT} ${DIFFCOEF} ${NX1}
else
    echo "Compilation error"
fi

python ../src/plotter.py

/opt/homebrew/Cellar/ffmpeg/5.1.2/bin/ffmpeg -r ${MOVIEFPS} -f image2 -i plot.%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p ${RUNNAME}.mp4

mv *.csv *.png outputs/.