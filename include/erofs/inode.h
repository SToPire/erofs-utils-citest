/* SPDX-License-Identifier: GPL-2.0+ OR Apache-2.0 */
/*
 * Copyright (C) 2018-2019 HUAWEI, Inc.
 *             http://www.huawei.com/
 * Created by Li Guifu <bluce.liguifu@huawei.com>
 * with heavy changes by Gao Xiang <xiang@kernel.org>
 */
#ifndef __EROFS_INODE_H
#define __EROFS_INODE_H

#ifdef __cplusplus
extern "C"
{
#endif

#include "erofs/internal.h"

static inline struct erofs_inode *erofs_igrab(struct erofs_inode *inode)
{
	++inode->i_count;
	return inode;
}

u32 erofs_new_encode_dev(dev_t dev);
unsigned char erofs_mode_to_ftype(umode_t mode);
umode_t erofs_ftype_to_mode(unsigned int ftype, unsigned int perm);
unsigned char erofs_ftype_to_dtype(unsigned int filetype);
void erofs_inode_manager_init(void);
void erofs_insert_ihash(struct erofs_inode *inode, dev_t dev, ino_t ino);
struct erofs_inode *erofs_iget(dev_t dev, ino_t ino);
struct erofs_inode *erofs_iget_by_nid(erofs_nid_t nid);
unsigned int erofs_iput(struct erofs_inode *inode);
erofs_nid_t erofs_lookupnid(struct erofs_inode *inode);
struct erofs_dentry *erofs_d_alloc(struct erofs_inode *parent,
				   const char *name);
int erofs_rebuild_dump_tree(struct erofs_inode *dir);
int erofs_init_empty_dir(struct erofs_inode *dir);
int __erofs_fill_inode(struct erofs_inode *inode, struct stat *st,
		       const char *path);
struct erofs_inode *erofs_new_inode(void);
struct erofs_inode *erofs_mkfs_build_tree_from_path(const char *path);
struct erofs_inode *erofs_mkfs_build_special_from_fd(int fd, const char *name);

#ifdef EROFS_MT_ENABLED
struct erofs_inode_fifo {
	pthread_mutex_t lock;
	pthread_cond_t full, empty;

	void *buf;

	size_t size, elem_size;
	size_t head, tail;
};

struct erofs_inode_fifo *erofs_alloc_inode_fifo(size_t size, size_t elem_size);
void erofs_push_inode_fifo(struct erofs_inode_fifo *q, void *elem);
void *erofs_pop_inode_fifo(struct erofs_inode_fifo *q);
void erofs_destroy_inode_fifo(struct erofs_inode_fifo *q);
#endif

#ifdef __cplusplus
}
#endif

#endif
