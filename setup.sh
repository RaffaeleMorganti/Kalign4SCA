# download kalign sources
wget https://github.com/TimoLassmann/kalign/archive/refs/tags/v3.4.0.tar.gz
tar -zxvf v3.4.0.tar.gz

# apply sca customization
mv custom/*.c kalign-3.4.0/lib/src/

# build & install
mkdir kalign-3.4.0/build && cd kalign-3.4.0/build
cmake ..
make
make install
cd ../..
