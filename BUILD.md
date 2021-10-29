
# Build instructions

## CMakeLists project

26.10.2021: for now, a few libs must be hard copied in `/lib`, we'll add them as submodule later. (AntTweakBar, field_tracing)

```
git submodule update --init --recursive

# Compile AntTweakBar
cd lib/AntTweakBar/src
make -j
cd ../../..

mkdir build
cd build
cmake ..
make -j parafashion
```

Then run with e.g. `./parafashion ../data/leggins/leggins.ply ../data/leggins/leggins.ply ../data/leggins/running_leggins_frames.txt`

Some other copy-pastable inputs: 
```
./parafashion ../data/gargo/gargoyle_rem.ply ../data/gargo/gargoyle_rem.ply
./parafashion ../data/wet/wet.obj ../data/wet/wet.obj

```

## From Nico's Parafashion

First, uncomment `LIBS += -lGLULIBS += -lGLU` in `Parafashion.pro`.

From the project folder:

```
sudo apt install qt5-default libqt5svg5-dev freeglut3-dev

# Compile AntTweakBar
cd lib/AntTweakBar1.16/src
make
cd ../../..

# temporary solution, set your own path here
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/corentin/Parafashion/lib/AntTweakBar/lib/
```

Then compile and run:
```
qmake Parafashion
make -j parafashion && ./parafashion test_data/leggins/leggins.ply test_data/leggins/leggins.ply 
```
