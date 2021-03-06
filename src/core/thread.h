/*
  Copyright (c) 2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef THREAD_H
#define THREAD_H

#include "core/error_api.h"

extern unsigned int gt_jobs; /* number of parallel threads to be used */

typedef struct GtThread GtThread;
typedef struct GtRWLock GtRWLock;
typedef struct GtMutex GtMutex;

typedef void* (*GtThreadFunc)(void *data);

/* Execute <function> (with <data> passed to it) in <gt_jobs> many parallel
   threads, if threading is enabled. Otherwise <function> is executed <gt_jobs>
   many times sequentially. */
int       gt_multithread(GtThreadFunc function, void *data, GtError *err);

/* Create a new thread which executes the given <function> (with <data> passed
   to it). Returns a <GtThread*> handle to the newly created thread, if
   successful. Returns NULL and sets <err> accordingly upon failure.  */
GtThread* gt_thread_new(GtThreadFunc function, void *data, GtError *err);

/* Delete the given <thread> handle. Does not stop the thread itself! */
void      gt_thread_delete(GtThread *thread);

void      gt_thread_join(GtThread *thread);

/* Return a new <GtRWLock*> object. */
GtRWLock* gt_rwlock_new(void);

/* Delete the given <rwlock>. */
void      gt_rwlock_delete(GtRWLock *rwlock);

/* Acquire a read lock for <rwlock>. */
#ifdef GT_THREADS_ENABLED
#define   gt_rwlock_rdlock(rwlock) \
          gt_rwlock_rdlock_func(rwlock)
void      gt_rwlock_rdlock_func(GtRWLock *rwlock);
#else
#define   gt_rwlock_rdlock(rwlock) \
          ((void) 0)
#endif

/* Acquire a write lock for <rwlock>. */
#ifdef GT_THREADS_ENABLED
#define   gt_rwlock_wrlock(rwlock) \
          gt_rwlock_wrlock_func(rwlock)
void      gt_rwlock_wrlock_func(GtRWLock *rwlock);
#else
#define   gt_rwlock_wrlock(rwlock) \
          ((void) 0)
#endif

/* Unlock the given <rwlock>. */
#ifdef GT_THREADS_ENABLED
#define   gt_rwlock_unlock(rwlock) \
          gt_rwlock_unlock_func(rwlock)
void      gt_rwlock_unlock_func(GtRWLock *rwlock);
#else
#define   gt_rwlock_unlock(rwlock) \
          ((void) 0)
#endif

/* Return a new <GtMutex*> object. */
GtMutex*  gt_mutex_new(void);

/* Delete the given <mutex>. */
void      gt_mutex_delete(GtMutex *mutex);

/* Lock the given <mutex>. */
#ifdef GT_THREADS_ENABLED
#define   gt_mutex_lock(mutex) \
          gt_mutex_lock_func(mutex)
void      gt_mutex_lock_func(GtMutex *mutex);
#else
#define   gt_mutex_lock(mutex) \
          ((void) 0)
#endif

/* Unlock the given <mutex>. */
#ifdef GT_THREADS_ENABLED
#define   gt_mutex_unlock(mutex) \
          gt_mutex_unlock_func(mutex)
void      gt_mutex_unlock_func(GtMutex *mutex);
#else
#define   gt_mutex_unlock(mutex) \
          ((void) 0)
#endif

#endif
