# parafashion

## BUILD

The following packages are required: `qt5-default libqt5svg5-dev freeglut3-dev` (e.g. `sudo apt install ...` on Debian based systems).

```
git clone https://github.com/nicopietroni/parafashion.git
git submodule update --init --recursive
mkdir build
cd build
cmake ..
make -j parafashion
```

Then, simply run parafashion with an obj or ply mesh as command-line argument: 
```
./parafashion ../data/leggins/leggins.ply
```


## (Alternatively) qmake compilation

From the project folder:
```
git submodule update --init --recursive
# Compile AntTweakBar
cd lib/AntTweakBar1.16/src
make
cd ../../..

#  set your own path here
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/[path_to_project]/lib/AntTweakBar/lib/
```

Then compile and run:
```
qmake parafashion
make -j parafashion 
./parafashion ../data/leggins/leggins.ply
```
