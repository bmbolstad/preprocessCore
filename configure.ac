dnl
dnl Configuration things for preprocessCore
dnl


AC_INIT("DESCRIPTION")

dnl
dnl Are things (still) the same ?
dnl (taken from the 'writing R extensions manual')

dnl Now find the compiler and compiler flags to use
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
       echo "could not determine R_HOME"
       exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CPP=`"${R_HOME}/bin/R" CMD config CPP`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
AC_PROG_CC
AC_PROG_CPP


AC_LANG(C)

# Checks for libraries.
use_pthreads=no
AC_SEARCH_LIBS([pthread_create], [pthread],
  [use_pthreads=yes])

# Checks for header files.
m4_warn([obsolete],
[The preprocessor macro 'STDC_HEADERS' is obsolete.
  Except in unusual embedded environments, you can safely include all
  C89 headers unconditionally.])dnl
# Autoupdate added the next two lines to ensure that your configure
# script's behavior did not change.  They are probably safe to remove.
AC_CHECK_INCLUDES_DEFAULT
AC_PROG_EGREP

AC_CHECK_HEADERS([stdlib.h])

if test "x${have_pthreads}" = xyes; then
  AC_CHECK_HEADERS([pthread.h], [],
    [use_pthreads=no])
fi


AC_MSG_CHECKING([if PTHREAD_STACK_MIN is defined])
AC_COMPILE_IFELSE([AC_LANG_SOURCE([
	#include <pthread.h>
	#include <limits.h>

      int main () {size_t stacksize = PTHREAD_STACK_MIN + 0x4000;
}     	  
    ])],[use_pthread_stack_min=yes], [use_pthread_stack_min=no])
AC_MSG_RESULT($use_pthread_stack_min)

if test "x$use_pthread_stack_min" = xno; then
	use_pthreads=no
fi


AC_MSG_CHECKING([if __pthread_get_minstack can be used])
AC_RUN_IFELSE([AC_LANG_SOURCE([
	#include <pthread.h>
	#include <limits.h> 
	#include <dlfcn.h>


	typedef size_t (*GetMinStack)(const pthread_attr_t *attr);

	GetMinStack _get_minstack_func = NULL;

	static void get_minstack_init() {
  		_get_minstack_func =
        (GetMinStack)dlsym(RTLD_DEFAULT, "__pthread_get_minstack");
	}

	
    int main () {
		get_minstack_init();
		if (_get_minstack_func == NULL){
			return 1;
		} else {
			return 0;
		}
	}

    ])],[use_pthread_get_minstack=yes], [use_pthread_get_minstack=no])
AC_MSG_RESULT($use_pthread_get_minstack)



dnl AC_MSG_CHECKING([if R is using flexiblas])
dnl BLAS_LIBS=`"${R_HOME}/bin/R" CMD config BLAS_LIBS`
dnl if echo "${BLAS_LIBS}" | grep flexiblas > /dev/null; then
dnl   AC_MSG_RESULT([flexiblas found. preprocessCore specific threading will be disabled])
dnl   use_pthreads=no
dnl else
dnl    AC_MSG_RESULT([flexiblas not found. preprocessCore threading will not be disabled])
dnl fi


dnl AC_MSG_CHECKING([if R is using openblasp])
dnl BLAS_LIBS=`"${R_HOME}/bin/R" CMD config BLAS_LIBS`
dnl if echo "${BLAS_LIBS}" | grep openblasp > /dev/null; then
dnl     AC_MSG_RESULT([openblasp found. preprocessCore specific threading will be disabled])
dnl     use_pthreads=no
dnl else
dnl     AC_MSG_RESULT([openblasp not found. preprocessCore threading will not be disabled])
dnl fi




AC_ARG_ENABLE([threading],
	AS_HELP_STRING([--disable-threading],[Disable threading]))



AS_IF([test "x$enable_threading" != "xno" ],[
	    if test "x${use_pthreads}" = "xno"; then
	       AC_MSG_NOTICE(------------------------------------------)
	       AC_MSG_NOTICE( Unable to find pthreads on this system.  )
	       AC_MSG_NOTICE( Building a single-threaded version.      )
	       AC_MSG_NOTICE(------------------------------------------)
	    fi


	    if test "x${use_pthreads}" = "xyes"; then
  	       AC_MSG_NOTICE(Enabling threading for preprocessCore)
  	       AC_DEFINE(USE_PTHREADS, 1)
		   if test "x${use_pthread_get_minstack}" = "xyes"; then
		    	AC_MSG_NOTICE(Using __pthread_get_minstack)
			    AC_DEFINE(INFER_MIN_STACKSIZE,1)
		   fi

	    fi
	    ],
	    [
	    AC_MSG_NOTICE(Disabling threading for preprocessCore)
	    ])





AC_OUTPUT(src/Makevars)
