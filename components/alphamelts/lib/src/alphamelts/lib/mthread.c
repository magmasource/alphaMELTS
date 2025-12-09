const char *mthread_ver(void) { return "$Id: mthread.c,v 1.1 2006/10/20 00:59:22 ghiorso Exp $"; }
/*
MELTS Source Code: RCS $Log: mthread.c,v $
MELTS Source Code: RCS Revision 1.1  2006/10/20 00:59:22  ghiorso
MELTS Source Code: RCS (1) Made initial modifications for thread safe code.
MELTS Source Code: RCS (2) Added support for XML I/O in batch mode
MELTS Source Code: RCS (3) Added support for Melts-batch listener for eventual integration into VIGMCS
MELTS Source Code: RCS
MELTS Source Code: RCS Revision 1.1  1999/11/13 22:09:00  ghiorso
MELTS Source Code: RCS Initial revision
MELTS Source Code: RCS
*/

/*
**++
**  FACILITY:  Thread Safe Silicate Melts Crystallization Package
**
**  MODULE DESCRIPTION:
**
**  Thread emulation routines
**  (file: MTHREAD.C)
**--
*/

#ifdef USE_PTHREADS

void MTHREAD_CLEANUPSPECIFIC(void) { ; }

#else

#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "mthread.h"

#ifdef DEBUG_THREAD
#undef DEBUG_THREAD
#endif 

#define	MTHREAD_DATAKEYS_MAX  512
#define OK                      0
#define MTHREAD_NEEDS_INIT      0
#define MTHREAD_DONE_INIT       1

int threadID;
char *threadList[MAX_CONNECTION_THREADS];

struct MTHREAD_KEY {
  long count;
  void (*destructor)(void *);
  void *specific_data[MAX_CONNECTION_THREADS];
};

static struct MTHREAD_KEY key_table[MTHREAD_DATAKEYS_MAX];

/* ==========================================================================
 * MTHREAD_KEY_CREATE()
 */
int MTHREAD_KEY_CREATE(MTHREAD_KEY_T *key, void (*destructor)(void *)) {
  for ((*key) = 0; (*key) < MTHREAD_DATAKEYS_MAX; (*key)++) {
    if (key_table[(*key)].count == 0) {
      key_table[(*key)].count++;
      key_table[(*key)].destructor = destructor;
#ifdef DEBUG_THREAD
      printf("Created key %d\n", *key);
#endif
      return(OK);
    }
  }
  return(EAGAIN);
}

/* ==========================================================================
 * MTHREAD_CLEANUPSPECIFIC()
 */
void MTHREAD_CLEANUPSPECIFIC(void) {
  void * data;
  int key;

  for (key = 0; key < MTHREAD_DATAKEYS_MAX; key++) {
    if (key_table[key].count) {
      if (key_table[key].specific_data[threadID]) {
#ifdef DEBUG_THREAD
        printf("Storage deleted for key %d in thread %d\n", key, threadID);
#endif
  	data = (void *) key_table[key].specific_data[threadID];
  	key_table[key].specific_data[threadID] = NULL;
  	if (key_table[key].destructor) key_table[key].destructor(data);
  	key_table[key].count--;
      }
    }
  }
}

/* ==========================================================================
 * MTHREAD_SETSPECIFIC()
 */
int MTHREAD_SETSPECIFIC(MTHREAD_KEY_T key, const void * value) {
  int ret;

#ifdef DEBUG_THREAD
  printf("Call to MTHREAD_SETSPECIFIC: key = %d, threadID = %d, value = %p\n", key, threadID, value);
#endif
  if (key < MTHREAD_DATAKEYS_MAX) {
    if (key_table[key].count) {
      if (key_table[key].specific_data[threadID] == NULL) {
  	if (value != NULL) key_table[key].count++;
      } else {
  	if (value == NULL) key_table[key].count--;
      }
      key_table[key].specific_data[threadID] = (void *) value;
      ret = OK;
    } else {
      ret = EINVAL;
    }
  } else {
    ret = EINVAL;
  }
  return(ret);
}

/* ==========================================================================
 * MTHREAD_GETSPECIFIC()
 */
void * MTHREAD_GETSPECIFIC(MTHREAD_KEY_T key) {
  void *ret;

  if ((key_table[key].specific_data[threadID]) && (key < MTHREAD_DATAKEYS_MAX)) {
    if (key_table[key].count) {
      ret = (void *) key_table[key].specific_data[threadID];
    } else {
      ret = NULL;
    }
  } else {
    ret = NULL;
  }
  return(ret);
}

/* ==========================================================================
 * MTHREAD_ONCE()
 */
int MTHREAD_ONCE(MTHREAD_ONCE_T *once_control, void (*init_routine)(void)) {
  if (*once_control == MTHREAD_NEEDS_INIT) {
    init_routine();
    *once_control = MTHREAD_DONE_INIT;
  }
  return(OK);
}

/* ==========================================================================
 * MTHREAD_MUTEX_LOCK()
 */
int MTHREAD_MUTEX_LOCK(MTHREAD_MUTEX_T *mutex) {
  return(OK); /* dummy routine */
}

/* ==========================================================================
 * MTHREAD_MUTEX_UNLOCK()
 */
int MTHREAD_MUTEX_UNLOCK(MTHREAD_MUTEX_T *mutex) {
  return(OK); /* dummy routine */
}

#endif /* USE_PTHREADS */

/* ==========================================================================
 * End of MTHREAD emulation routines
 */
