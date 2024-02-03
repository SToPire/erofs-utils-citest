// SPDX-License-Identifier: GPL-2.0+
#include "erofs/err.h"
#include <stdlib.h>
#include "erofs/queue.h"

struct erofs_queue *erofs_queue_create(size_t size, size_t elem_size)
{
	struct erofs_queue *q = malloc(sizeof(*q));

	pthread_mutex_init(&q->lock, NULL);
	pthread_cond_init(&q->empty, NULL);
	pthread_cond_init(&q->full, NULL);

	q->size = size;
	q->elem_size = elem_size;
	q->head = 0;
	q->tail = 0;
	q->buf = calloc(size, elem_size);
	if (!q->buf)
		return ERR_PTR(-ENOMEM);

	return q;
}

void erofs_queue_push(struct erofs_queue *q, void *elem)
{
	pthread_mutex_lock(&q->lock);

	while ((q->tail + 1) % q->size == q->head)
		pthread_cond_wait(&q->full, &q->lock);

	memcpy(q->buf + q->tail * q->elem_size, elem, q->elem_size);
	q->tail = (q->tail + 1) % q->size;

	pthread_cond_signal(&q->empty);
	pthread_mutex_unlock(&q->lock);
}

void *erofs_queue_pop(struct erofs_queue *q)
{
    void *elem;

    pthread_mutex_lock(&q->lock);

    while (q->head == q->tail)
        pthread_cond_wait(&q->empty, &q->lock);

    elem = q->buf + q->head * q->elem_size;
    q->head = (q->head + 1) % q->size;

    pthread_cond_signal(&q->full);
    pthread_mutex_unlock(&q->lock);

    return elem;
}

void erofs_queue_destroy(struct erofs_queue *q)
{
	pthread_mutex_destroy(&q->lock);
	pthread_cond_destroy(&q->empty);
	pthread_cond_destroy(&q->full);
	free(q->buf);
	free(q);
}