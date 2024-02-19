/* SPDX-License-Identifier: GPL-2.0+ OR Apache-2.0 */
#ifndef __EROFS_WORKQUEUE_H
#define __EROFS_WORKQUEUE_H

#include "internal.h"

struct erofs_work;

typedef void erofs_wq_func_t(struct erofs_work *);
typedef void erofs_wq_priv_fini_t(void *);

struct erofs_work {
	void (*func)(struct erofs_work *work);
	struct erofs_work *next;
	void *priv;
};

struct erofs_workqueue {
	struct erofs_work *head, *tail;
	pthread_mutex_t lock;
	pthread_cond_t cond_empty;
	pthread_cond_t cond_full;
	pthread_t *workers;
	unsigned int nworker;
	unsigned int max_jobs;
	unsigned int job_count;
	bool shutdown;
	size_t priv_size;
	erofs_wq_priv_fini_t *priv_fini;
};

int erofs_alloc_workqueue(struct erofs_workqueue *wq, unsigned int nworker,
			 unsigned int max_jobs, size_t priv_size,
			 erofs_wq_priv_fini_t *priv_fini);
int erofs_queue_work(struct erofs_workqueue *wq, struct erofs_work *work);
int erofs_destroy_workqueue(struct erofs_workqueue *wq);
#endif