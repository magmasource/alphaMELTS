#ifndef _Mthread_h
#define _Mthread_h

/*
MELTS Source Code: RCS $Log: mthread.h,v $
MELTS Source Code: RCS Revision 1.1  2006/10/20 00:59:22  ghiorso
MELTS Source Code: RCS (1) Made initial modifications for thread safe code.
MELTS Source Code: RCS (2) Added support for XML I/O in batch mode
MELTS Source Code: RCS (3) Added support for Melts-batch listener for eventual integration into VIGMCS
MELTS Source Code: RCS
*/

/*
**++
**  FACILITY:  Thread Safe Silicate Melts Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      pThread emulation include file (file: MTHREAD.H)
**--
*/


#ifdef USE_PTHREADS

#include <pthread.h>

#define MTHREAD_KEY_T             pthread_key_t
#define MTHREAD_ONCE_T            pthread_once_t
#define MTHREAD_ONCE_INIT         PTHREAD_ONCE_INIT
#define MTHREAD_MUTEX_T           pthread_mutex_t
#define MTHREAD_MUTEX_INITIALIZER PTHREAD_MUTEX_INITIALIZER

#define MTHREAD_GETSPECIFIC(x)    pthread_getspecific((x))
#define MTHREAD_KEY_CREATE(x, y)  pthread_key_create((x), (y))
#define MTHREAD_ONCE(x, y)        pthread_once((x), (y))
#define MTHREAD_SETSPECIFIC(x, y) pthread_setspecific((x), (y))
#define MTHREAD_MUTEX_LOCK(x)     pthread_mutex_lock((x))
#define MTHREAD_MUTEX_UNLOCK(x)   pthread_mutex_unlock((x))

void  MTHREAD_CLEANUPSPECIFIC(void);

#else

typedef unsigned int MTHREAD_KEY_T;
typedef int          MTHREAD_ONCE_T;
typedef int          MTHREAD_MUTEX_T;

#define MTHREAD_ONCE_INIT         0
#define MTHREAD_MUTEX_INITIALIZER 0

void  MTHREAD_CLEANUPSPECIFIC(void); 
void *MTHREAD_GETSPECIFIC(MTHREAD_KEY_T key);
int   MTHREAD_KEY_CREATE(MTHREAD_KEY_T *key, void (*destructor)(void *));
int   MTHREAD_ONCE(MTHREAD_ONCE_T *once_control, void (*init_routine)(void));
int   MTHREAD_SETSPECIFIC(MTHREAD_KEY_T key, const void *value);
int   MTHREAD_MUTEX_LOCK(MTHREAD_MUTEX_T *mutex);
int   MTHREAD_MUTEX_UNLOCK(MTHREAD_MUTEX_T *mutex);

#define MAX_CONNECTION_THREADS 25

extern int threadID;
extern char *threadList[MAX_CONNECTION_THREADS];

#endif

#endif /* _Mthread_h */
