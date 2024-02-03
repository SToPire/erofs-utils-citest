/* SPDX-License-Identifier: GPL-2.0+ */
#ifndef __EROFS_QUEUE_H
#define __EROFS_QUEUE_H

#include "internal.h"

struct erofs_queue {
    pthread_mutex_t lock;
    pthread_cond_t full, empty;

    void *buf;

    size_t size, elem_size;
    size_t head, tail;
};

struct erofs_queue* erofs_queue_create(size_t size, size_t elem_size);
void erofs_queue_push(struct erofs_queue *q, void *elem);
void *erofs_queue_pop(struct erofs_queue *q);
void erofs_queue_destroy(struct erofs_queue *q);

#endif