# neccessary PERL dependencies to install MEME
# cpan -u

apt-get remove -y perl-homedir
cpan install File::Which HTML::Template JSON XML::Parser Sys::Info
