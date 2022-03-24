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