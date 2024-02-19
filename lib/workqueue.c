// SPDX-License-Identifier: GPL-2.0+ OR Apache-2.0
#include <pthread.h>
#include <stdlib.h>
#include "erofs/workqueue.h"

static void *worker_thread(void *arg)
{
	struct erofs_workqueue *wq = arg;
	struct erofs_work *work;
	void *priv = NULL;

	if (wq->on_start)
		priv = (wq->on_start)(wq, NULL);

	while (true) {
		pthread_mutex_lock(&wq->lock);

		while (wq->job_count == 0 && !wq->shutdown)
			pthread_cond_wait(&wq->cond_empty, &wq->lock);
		if (wq->job_count == 0 && wq->shutdown) {
			pthread_mutex_unlock(&wq->lock);
			break;
		}

		work = wq->head;
		wq->head = work->next;
		if (!wq->head)
			wq->tail = NULL;
		wq->job_count--;

		if (wq->job_count == wq->max_jobs - 1)
			pthread_cond_broadcast(&wq->cond_full);

		pthread_mutex_unlock(&wq->lock);

		work->priv = priv;
		work->func(work);
	}

	if (wq->on_exit)
		(wq->on_exit)(wq, priv);

	return NULL;
}

int erofs_alloc_workqueue(struct erofs_workqueue *wq, unsigned int nworker,
			  unsigned int max_jobs, erofs_wq_func_t start_func,
			  erofs_wq_func_t exit_func)
{
	unsigned int i;
	int ret;

	if (!wq || nworker <= 0 || max_jobs <= 0)
		return -EINVAL;

	wq->head = wq->tail = NULL;
	wq->nworker = nworker;
	wq->max_jobs = max_jobs;
	wq->job_count = 0;
	wq->shutdown = false;
	wq->on_start = start_func;
	wq->on_exit = exit_func;
	pthread_mutex_init(&wq->lock, NULL);
	pthread_cond_init(&wq->cond_empty, NULL);
	pthread_cond_init(&wq->cond_full, NULL);

	wq->workers = malloc(nworker * sizeof(pthread_t));
	if (!wq->workers)
		return -ENOMEM;

	for (i = 0; i < nworker; i++) {
		ret = pthread_create(&wq->workers[i], NULL, worker_thread, wq);
		if (ret) {
			while (i)
				pthread_cancel(wq->workers[--i]);
			free(wq->workers);
			return ret;
		}
	}

	return 0;
}

int erofs_queue_work(struct erofs_workqueue *wq, struct erofs_work *work)
{
	if (!wq || !work)
		return -EINVAL;

	pthread_mutex_lock(&wq->lock);

	while (wq->job_count == wq->max_jobs)
		pthread_cond_wait(&wq->cond_full, &wq->lock);

	work->next = NULL;
	if (!wq->head)
		wq->head = work;
	else
		wq->tail->next = work;
	wq->tail = work;
	wq->job_count++;

	pthread_cond_signal(&wq->cond_empty);
	pthread_mutex_unlock(&wq->lock);

	return 0;
}

int erofs_destroy_workqueue(struct erofs_workqueue *wq)
{
	unsigned int i;

	if (!wq)
		return -EINVAL;

	pthread_mutex_lock(&wq->lock);
	wq->shutdown = true;
	pthread_cond_broadcast(&wq->cond_empty);
	pthread_mutex_unlock(&wq->lock);

	for (i = 0; i < wq->nworker; i++)
		pthread_join(wq->workers[i], NULL);

	free(wq->workers);
	pthread_mutex_destroy(&wq->lock);
	pthread_cond_destroy(&wq->cond_empty);
	pthread_cond_destroy(&wq->cond_full);

	return 0;
}