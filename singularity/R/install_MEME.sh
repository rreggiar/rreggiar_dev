# pulled from the website
# downloads and installs MEME-suite to opt/meme-5.5.1
# need to set this path in R; it is not default for {memes::}

version=5.5.1
wget http://meme-suite.org/meme-software/$version/meme-$version.tar.gz
tar zxf meme-$version.tar.gz
cd meme-$version
./configure --prefix=/etc/meme-5.5.1
make
make test
make install
