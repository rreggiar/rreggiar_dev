# neccessary PERL dependencies to install MEME
# cpan -u

apt-get remove -y perl-homedir
cpan install File::Which HTML::Template JSON XML::Parser Sys::Info Log::Log4perl Math::CDF XML::Compile::SOAP11 XML::Compile::WSDL11
