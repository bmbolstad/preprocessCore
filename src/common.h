#ifndef COMMON_H
#define COMMON_H 1

#ifdef USE_PTHREADS
#include <pthread.h>
extern pthread_mutex_t mutex_R;
#endif

#endif
