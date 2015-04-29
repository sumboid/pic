AC_INIT([poisson],[0.1],[ilya.sumb@gmail.com])
AM_INIT_AUTOMAKE([foreign silent-rules])
AM_SILENT_RULES([yes])

AC_PROG_CXX([mpic++])

CXXFLAGS="-O0 -g -std=c++11 -I/home/yoba/.local/ts/include -Wall -Wextra -Werror"
LDFLAGS="-L/home/yoba/.local/ts/lib"
LIBS="-lts"

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
