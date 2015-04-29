#!/usr/bin/env sh

list="aclocal.m4 \
			autom4te.cache \
			config.guess \
			config.log \
		  config.status \
		  config.sub \
		  configure \
		  depcomp \
		  install-sh \
		  .libs \
		  libtool \
		  Makefile \
		  Makefile.in \
		  missing \
		  m4/libtool.m4 \
		  m4/lt~obsolete.m4 \
		  m4/ltoptions.m4 \
		  m4/ltsugar.m4 \
		  m4/ltversion.m4"	

generate () {
	aclocal -I m4 --install &&
	autoconf --force &&
	libtoolize &&
	automake --add-missing 
}

clean () {
	rm -rf $list
}

clean
generate
