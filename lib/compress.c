// SPDX-License-Identifier: GPL-2.0+ OR Apache-2.0
/*
 * Copyright (C) 2018-2019 HUAWEI, Inc.
 *             http://www.huawei.com/
 * Created by Miao Xie <miaoxie@huawei.com>
 * with heavy changes by Gao Xiang <gaoxiang25@huawei.com>
 */
#ifndef _LARGEFILE64_SOURCE
#define _LARGEFILE64_SOURCE
#endif
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include "erofs/print.h"
#include "erofs/io.h"
#include "erofs/cache.h"
#include "erofs/compress.h"
#include "erofs/dedupe.h"
#include "compressor.h"
#include "erofs/block_list.h"
#include "erofs/compress_hints.h"
#include "erofs/fragments.h"
#ifdef EROFS_MT_ENABLED
#include "erofs/workqueue.h"
#endif
#ifdef HAVE_LINUX_FALLOC_H
#include <linux/falloc.h>
#endif

#if defined(HAVE_FALLOCATE) && defined(FALLOC_FL_PUNCH_HOLE)
#define USE_PER_WORKER_TMPFILE 1
#endif

/* compressing configuration specified by users */
struct erofs_compress_cfg {
	struct erofs_compress handle;
	unsigned int algorithmtype;
	bool enable;
} erofs_ccfg[EROFS_MAX_COMPR_CFGS];

struct z_erofs_extent_item {
	struct list_head list;
	struct z_erofs_inmem_extent e;
};

struct z_erofs_file_compress_ctx {
	struct erofs_inode *inode;
	int fd;
	unsigned int pclustersize;

	u32 tof_chksum;
	bool fix_dedupedfrag;
	bool fragemitted;
};

struct z_erofs_vle_compress_ctx {
	struct z_erofs_file_compress_ctx *fctx;

	u8 *queue;
	struct list_head extents;
	struct z_erofs_extent_item *pivot;

	struct erofs_compress *chandle;
	char *destbuf;

	unsigned int head, tail;
	erofs_off_t remaining;
	erofs_blk_t blkaddr;		/* pointing to the next blkaddr */
	erofs_blk_t compressed_blocks;
	u16 clusterofs;

	int seg_num, seg_idx;
	FILE *tmpfile;
	off_t tmpfile_off;
};

struct z_erofs_write_index_ctx {
	struct erofs_inode *inode;
	struct list_head *extents;
	u16 clusterofs;
	erofs_blk_t blkaddr, blkoff;
	u8 *metacur;
};

#ifdef EROFS_MT_ENABLED
struct erofs_compress_wq_private {
	bool init;
	u8 *queue;
	char *destbuf;
	struct erofs_compress_cfg *ccfg;
	FILE* tmpfile;
};

struct erofs_compress_work {
	/* Note: struct erofs_work must be the first member */
	struct erofs_work work;
	struct z_erofs_vle_compress_ctx ctx;

	unsigned int alg_id;
	char *alg_name;
	unsigned int comp_level;
	unsigned int dict_size;

	int ret;

	struct erofs_compress_work *next;
};

static struct {
	struct erofs_workqueue wq;
	struct erofs_compress_work *idle;
	pthread_mutex_t mutex;
	pthread_cond_t cond;
	int nfini;
} z_erofs_mt_ctrl;
#endif

static bool z_erofs_mt_enabled;
static u8 *z_erofs_global_queue;

#define Z_EROFS_LEGACY_MAP_HEADER_SIZE	Z_EROFS_FULL_INDEX_ALIGN(0)

static void z_erofs_write_indexes_final(struct z_erofs_write_index_ctx *ctx)
{
	const unsigned int type = Z_EROFS_LCLUSTER_TYPE_PLAIN;
	struct z_erofs_lcluster_index di;

	if (!ctx->clusterofs)
		return;

	di.di_clusterofs = cpu_to_le16(ctx->clusterofs);
	di.di_u.blkaddr = 0;
	di.di_advise = cpu_to_le16(type << Z_EROFS_LI_LCLUSTER_TYPE_BIT);

	memcpy(ctx->metacur, &di, sizeof(di));
	ctx->metacur += sizeof(di);
}

static void z_erofs_write_extent(struct z_erofs_write_index_ctx *ctx,
				 struct z_erofs_inmem_extent *e)
{
	struct erofs_inode *inode = ctx->inode;
	struct erofs_sb_info *sbi = inode->sbi;
	unsigned int clusterofs = ctx->clusterofs;
	unsigned int count = e->length;
	unsigned int d0 = 0, d1 = (clusterofs + count) / erofs_blksiz(sbi);
	struct z_erofs_lcluster_index di;
	unsigned int type, advise;

	DBG_BUGON(!count);
	di.di_clusterofs = cpu_to_le16(ctx->clusterofs);

	/* whether the tail-end (un)compressed block or not */
	if (!d1) {
		/*
		 * A lcluster cannot have three parts with the middle one which
		 * is well-compressed for !ztailpacking cases.
		 */
		DBG_BUGON(!e->raw && !cfg.c_ztailpacking && !cfg.c_fragments);
		DBG_BUGON(e->partial);
		type = e->raw ? Z_EROFS_LCLUSTER_TYPE_PLAIN :
			Z_EROFS_LCLUSTER_TYPE_HEAD1;
		advise = type << Z_EROFS_LI_LCLUSTER_TYPE_BIT;
		di.di_advise = cpu_to_le16(advise);

		if (inode->datalayout == EROFS_INODE_COMPRESSED_FULL &&
		    !e->compressedblks) {
			di.di_u.blkaddr = cpu_to_le32(inode->fragmentoff >> 32);
		} else if (z_erofs_mt_enabled) {
			di.di_u.blkaddr =
				cpu_to_le32(ctx->blkaddr + ctx->blkoff);
			ctx->blkoff += e->compressedblks;
		} else {
			di.di_u.blkaddr = cpu_to_le32(e->blkaddr);
		}
		memcpy(ctx->metacur, &di, sizeof(di));
		ctx->metacur += sizeof(di);

		/* don't add the final index if the tail-end block exists */
		ctx->clusterofs = 0;
		return;
	}

	do {
		advise = 0;
		/* XXX: big pcluster feature should be per-inode */
		if (d0 == 1 && erofs_sb_has_big_pcluster(sbi)) {
			type = Z_EROFS_LCLUSTER_TYPE_NONHEAD;
			di.di_u.delta[0] = cpu_to_le16(e->compressedblks |
						       Z_EROFS_LI_D0_CBLKCNT);
			di.di_u.delta[1] = cpu_to_le16(d1);
		} else if (d0) {
			type = Z_EROFS_LCLUSTER_TYPE_NONHEAD;

			/*
			 * If the |Z_EROFS_VLE_DI_D0_CBLKCNT| bit is set, parser
			 * will interpret |delta[0]| as size of pcluster, rather
			 * than distance to last head cluster. Normally this
			 * isn't a problem, because uncompressed extent size are
			 * below Z_EROFS_VLE_DI_D0_CBLKCNT * BLOCK_SIZE = 8MB.
			 * But with large pcluster it's possible to go over this
			 * number, resulting in corrupted compressed indices.
			 * To solve this, we replace d0 with
			 * Z_EROFS_VLE_DI_D0_CBLKCNT-1.
			 */
			if (d0 >= Z_EROFS_LI_D0_CBLKCNT)
				di.di_u.delta[0] = cpu_to_le16(
						Z_EROFS_LI_D0_CBLKCNT - 1);
			else
				di.di_u.delta[0] = cpu_to_le16(d0);
			di.di_u.delta[1] = cpu_to_le16(d1);
		} else {
			type = e->raw ? Z_EROFS_LCLUSTER_TYPE_PLAIN :
				Z_EROFS_LCLUSTER_TYPE_HEAD1;

			if (inode->datalayout == EROFS_INODE_COMPRESSED_FULL &&
			    !e->compressedblks) {
				di.di_u.blkaddr = cpu_to_le32(inode->fragmentoff >> 32);
			} else if (z_erofs_mt_enabled) {
				di.di_u.blkaddr =
					cpu_to_le32(ctx->blkaddr + ctx->blkoff);
				ctx->blkoff += e->compressedblks;
			} else {
				di.di_u.blkaddr = cpu_to_le32(e->blkaddr);
			}

			if (e->partial) {
				DBG_BUGON(e->raw);
				advise |= Z_EROFS_LI_PARTIAL_REF;
			}
		}
		advise |= type << Z_EROFS_LI_LCLUSTER_TYPE_BIT;
		di.di_advise = cpu_to_le16(advise);

		memcpy(ctx->metacur, &di, sizeof(di));
		ctx->metacur += sizeof(di);

		count -= erofs_blksiz(sbi) - clusterofs;
		clusterofs = 0;

		++d0;
		--d1;
	} while (clusterofs + count >= erofs_blksiz(sbi));

	ctx->clusterofs = clusterofs + count;
}

static void z_erofs_write_indexes(struct z_erofs_write_index_ctx *ctx)
{
	struct z_erofs_extent_item *ei, *n;

	ctx->clusterofs = 0;
	list_for_each_entry_safe(ei, n, ctx->extents, list) {
		z_erofs_write_extent(ctx, &ei->e);

		list_del(&ei->list);
		free(ei);
	}
	z_erofs_write_indexes_final(ctx);
}

static bool z_erofs_need_refill(struct z_erofs_vle_compress_ctx *ctx)
{
	const bool final = !ctx->remaining;
	unsigned int qh_aligned, qh_after;
	struct erofs_inode *inode = ctx->fctx->inode;

	if (final || ctx->head < EROFS_CONFIG_COMPR_MAX_SZ)
		return false;

	qh_aligned = round_down(ctx->head, erofs_blksiz(inode->sbi));
	qh_after = ctx->head - qh_aligned;
	memmove(ctx->queue, ctx->queue + qh_aligned, ctx->tail - qh_aligned);
	ctx->tail -= qh_aligned;
	ctx->head = qh_after;
	return true;
}

static struct z_erofs_extent_item dummy_pivot = {
	.e.length = 0
};

static void z_erofs_commit_extent(struct z_erofs_vle_compress_ctx *ctx,
				  struct z_erofs_extent_item *ei)
{
	if (ei == &dummy_pivot)
		return;

	list_add_tail(&ei->list, &ctx->extents);
	ctx->clusterofs = (ctx->clusterofs + ei->e.length) &
			  (erofs_blksiz(ctx->fctx->inode->sbi) - 1);
}

static int z_erofs_compress_dedupe(struct z_erofs_vle_compress_ctx *ctx,
				   unsigned int *len)
{
	struct erofs_inode *inode = ctx->fctx->inode;
	const unsigned int lclustermask = (1 << inode->z_logical_clusterbits) - 1;
	struct erofs_sb_info *sbi = inode->sbi;
	struct z_erofs_extent_item *ei = ctx->pivot;

	if (!ei)
		return 0;

	/*
	 * No need dedupe for packed inode since it is composed of
	 * fragments which have already been deduplicated.
	 */
	if (erofs_is_packed_inode(inode))
		goto out;

	do {
		struct z_erofs_dedupe_ctx dctx = {
			.start = ctx->queue + ctx->head - ({ int rc;
				if (ei->e.length <= erofs_blksiz(sbi))
					rc = 0;
				else if (ei->e.length - erofs_blksiz(sbi) >= ctx->head)
					rc = ctx->head;
				else
					rc = ei->e.length - erofs_blksiz(sbi);
				rc; }),
			.end = ctx->queue + ctx->head + *len,
			.cur = ctx->queue + ctx->head,
		};
		int delta;

		if (z_erofs_dedupe_match(&dctx))
			break;

		delta = ctx->queue + ctx->head - dctx.cur;
		/*
		 * For big pcluster dedupe, leave two indices at least to store
		 * CBLKCNT as the first step.  Even laterly, an one-block
		 * decompresssion could be done as another try in practice.
		 */
		if (dctx.e.compressedblks > 1 &&
		    ((ctx->clusterofs + ei->e.length - delta) & lclustermask) +
			dctx.e.length < 2 * (lclustermask + 1))
			break;

		ctx->pivot = malloc(sizeof(struct z_erofs_extent_item));
		if (!ctx->pivot) {
			z_erofs_commit_extent(ctx, ei);
			return -ENOMEM;
		}

		if (delta) {
			DBG_BUGON(delta < 0);
			DBG_BUGON(!ei->e.length);

			/*
			 * For big pcluster dedupe, if we decide to shorten the
			 * previous big pcluster, make sure that the previous
			 * CBLKCNT is still kept.
			 */
			if (ei->e.compressedblks > 1 &&
			    (ctx->clusterofs & lclustermask) + ei->e.length
				- delta < 2 * (lclustermask + 1))
				break;
			ei->e.partial = true;
			ei->e.length -= delta;
		}

		/* fall back to noncompact indexes for deduplication */
		inode->z_advise &= ~Z_EROFS_ADVISE_COMPACTED_2B;
		inode->datalayout = EROFS_INODE_COMPRESSED_FULL;
		erofs_sb_set_dedupe(sbi);

		sbi->saved_by_deduplication +=
			dctx.e.compressedblks * erofs_blksiz(sbi);
		erofs_dbg("Dedupe %u %scompressed data (delta %d) to %u of %u blocks",
			  dctx.e.length, dctx.e.raw ? "un" : "",
			  delta, dctx.e.blkaddr, dctx.e.compressedblks);

		z_erofs_commit_extent(ctx, ei);
		ei = ctx->pivot;
		init_list_head(&ei->list);
		ei->e = dctx.e;

		ctx->head += dctx.e.length - delta;
		DBG_BUGON(*len < dctx.e.length - delta);
		*len -= dctx.e.length - delta;

		if (z_erofs_need_refill(ctx))
			return 1;
	} while (*len);
out:
	z_erofs_commit_extent(ctx, ei);
	ctx->pivot = NULL;
	return 0;
}

static int write_uncompressed_extent(struct z_erofs_vle_compress_ctx *ctx,
				     unsigned int len, char *dst)
{
	struct erofs_inode *inode = ctx->fctx->inode;
	struct erofs_sb_info *sbi = inode->sbi;
	unsigned int count = min(erofs_blksiz(sbi), len);
	unsigned int interlaced_offset, rightpart;
	int ret;

	/* write interlaced uncompressed data if needed */
	if (inode->z_advise & Z_EROFS_ADVISE_INTERLACED_PCLUSTER)
		interlaced_offset = ctx->clusterofs;
	else
		interlaced_offset = 0;
	rightpart = min(erofs_blksiz(sbi) - interlaced_offset, count);

	memset(dst, 0, erofs_blksiz(sbi));

	memcpy(dst + interlaced_offset, ctx->queue + ctx->head, rightpart);
	memcpy(dst, ctx->queue + ctx->head + rightpart, count - rightpart);

	if (ctx->tmpfile) {
		erofs_dbg("Writing %u uncompressed data to tmpfile", count);
		ret = fwrite(dst, erofs_blksiz(sbi), 1, ctx->tmpfile);
		if (ret != 1)
			return -EIO;
		fflush(ctx->tmpfile);
	} else {
		erofs_dbg("Writing %u uncompressed data to block %u", count,
			  ctx->blkaddr);
		ret = blk_write(sbi, dst, ctx->blkaddr, 1);
		if (ret)
			return ret;
	}
	return count;
}

static unsigned int z_erofs_get_max_pclustersize(struct erofs_inode *inode)
{
	unsigned int pclusterblks;

	if (erofs_is_packed_inode(inode))
		pclusterblks = cfg.c_pclusterblks_packed;
#ifndef NDEBUG
	else if (cfg.c_random_pclusterblks)
		pclusterblks = 1 + rand() % cfg.c_pclusterblks_max;
#endif
	else if (cfg.c_compress_hints_file) {
		z_erofs_apply_compress_hints(inode);
		DBG_BUGON(!inode->z_physical_clusterblks);
		pclusterblks = inode->z_physical_clusterblks;
	} else {
		pclusterblks = cfg.c_pclusterblks_def;
	}
	return pclusterblks * erofs_blksiz(inode->sbi);
}

static int z_erofs_fill_inline_data(struct erofs_inode *inode, void *data,
				    unsigned int len, bool raw)
{
	inode->z_advise |= Z_EROFS_ADVISE_INLINE_PCLUSTER;
	inode->idata_size = len;
	inode->compressed_idata = !raw;

	inode->idata = malloc(inode->idata_size);
	if (!inode->idata)
		return -ENOMEM;
	erofs_dbg("Recording %u %scompressed inline data",
		  inode->idata_size, raw ? "un" : "");
	memcpy(inode->idata, data, inode->idata_size);
	return len;
}

static void tryrecompress_trailing(struct z_erofs_vle_compress_ctx *ctx,
				   struct erofs_compress *ec,
				   void *in, unsigned int *insize,
				   void *out, unsigned int *compressedsize)
{
	struct erofs_sb_info *sbi = ctx->fctx->inode->sbi;
	char tmp[Z_EROFS_PCLUSTER_MAX_SIZE];
	unsigned int count;
	int ret = *compressedsize;

	/* no need to recompress */
	if (!(ret & (erofs_blksiz(sbi) - 1)))
		return;

	count = *insize;
	ret = erofs_compress_destsize(ec, in, &count, (void *)tmp,
				      rounddown(ret, erofs_blksiz(sbi)));
	if (ret <= 0 || ret + (*insize - count) >=
			roundup(*compressedsize, erofs_blksiz(sbi)))
		return;

	/* replace the original compressed data if any gain */
	memcpy(out, tmp, ret);
	*insize = count;
	*compressedsize = ret;
}

static bool z_erofs_fixup_deduped_fragment(struct z_erofs_vle_compress_ctx *ctx,
					   unsigned int len)
{
	struct erofs_inode *inode = ctx->fctx->inode;
	struct erofs_sb_info *sbi = inode->sbi;
	const unsigned int newsize = ctx->remaining + len;

	DBG_BUGON(!inode->fragment_size);

	/* try to fix again if it gets larger (should be rare) */
	if (inode->fragment_size < newsize) {
		ctx->fctx->pclustersize =
			min_t(erofs_off_t, z_erofs_get_max_pclustersize(inode),
			      roundup(newsize - inode->fragment_size,
				      erofs_blksiz(sbi)));
		return false;
	}

	inode->fragmentoff += inode->fragment_size - newsize;
	inode->fragment_size = newsize;

	erofs_dbg("Reducing fragment size to %u at %llu",
		  inode->fragment_size, inode->fragmentoff | 0ULL);

	/* it's the end */
	DBG_BUGON(ctx->tail - ctx->head + ctx->remaining != newsize);
	ctx->head = ctx->tail;
	ctx->remaining = 0;
	return true;
}

static int __z_erofs_compress_one(struct z_erofs_vle_compress_ctx *ctx,
				  struct z_erofs_inmem_extent *e)
{
	static char
		global_dstbuf[EROFS_CONFIG_COMPR_MAX_SZ + EROFS_MAX_BLOCK_SIZE];
	char *dstbuf = ctx->destbuf ? ctx->destbuf : global_dstbuf;
	struct z_erofs_file_compress_ctx *fctx = ctx->fctx;
	struct erofs_inode *inode = fctx->inode;
	struct erofs_sb_info *sbi = inode->sbi;
	unsigned int blksz = erofs_blksiz(sbi);
	char *const dst = dstbuf + blksz;
	struct erofs_compress *const h = ctx->chandle;
	unsigned int len = ctx->tail - ctx->head;
	bool is_packed_inode = erofs_is_packed_inode(inode);
	bool final = !ctx->remaining;
	bool may_packing = (cfg.c_fragments && final && !is_packed_inode &&
			    !z_erofs_mt_enabled);
	bool may_inline = (cfg.c_ztailpacking && final && !may_packing &&
			   !z_erofs_mt_enabled);
	unsigned int compressedsize;
	int ret;

	if (len <= fctx->pclustersize) {
		if (!final || !len)
			return 1;
		if (may_packing) {
			if (inode->fragment_size && !fctx->fix_dedupedfrag) {
				fctx->pclustersize = roundup(len, blksz);
				goto fix_dedupedfrag;
			}
			e->length = len;
			goto frag_packing;
		}
		if (!may_inline && len <= blksz)
			goto nocompression;
	}

	e->length = min(len, cfg.c_max_decompressed_extent_bytes);
	ret = erofs_compress_destsize(h, ctx->queue + ctx->head,
				      &e->length, dst, fctx->pclustersize);
	if (ret <= 0) {
		erofs_err("failed to compress %s: %s", inode->i_srcpath,
			  erofs_strerror(ret));
		return ret;
	}

	compressedsize = ret;
	/* even compressed size is smaller, there is no real gain */
	if (!(may_inline && e->length == len && ret < blksz))
		ret = roundup(ret, blksz);

	/* check if there is enough gain to keep the compressed data */
	if (ret * h->compress_threshold / 100 >= e->length) {
		if (may_inline && len < blksz) {
			ret = z_erofs_fill_inline_data(inode,
					ctx->queue + ctx->head, len, true);
		} else {
			may_inline = false;
			may_packing = false;
nocompression:
			/* TODO: reset clusterofs to 0 if permitted */
			ret = write_uncompressed_extent(ctx, len, dst);
		}

		if (ret < 0)
			return ret;
		e->length = ret;

		/*
		 * XXX: For now, we have to leave `ctx->compressedblk = 1'
		 * since there is no way to generate compressed indexes after
		 * the time that ztailpacking is decided.
		 */
		e->compressedblks = 1;
		e->raw = true;
	} else if (may_packing && len == e->length &&
		   compressedsize < fctx->pclustersize &&
		   (!inode->fragment_size || fctx->fix_dedupedfrag)) {
frag_packing:
		ret = z_erofs_pack_fragments(inode, ctx->queue + ctx->head,
					     len, fctx->tof_chksum);
		if (ret < 0)
			return ret;
		e->compressedblks = 0; /* indicate a fragment */
		e->raw = false;
		fctx->fragemitted = true;
	/* tailpcluster should be less than 1 block */
	} else if (may_inline && len == e->length && compressedsize < blksz) {
		if (ctx->clusterofs + len <= blksz) {
			inode->eof_tailraw = malloc(len);
			if (!inode->eof_tailraw)
				return -ENOMEM;

			memcpy(inode->eof_tailraw, ctx->queue + ctx->head, len);
			inode->eof_tailrawsize = len;
		}

		ret = z_erofs_fill_inline_data(inode, dst,
				compressedsize, false);
		if (ret < 0)
			return ret;
		e->compressedblks = 1;
		e->raw = false;
	} else {
		unsigned int tailused, padding;

		/*
		 * If there's space left for the last round when deduping
		 * fragments, try to read the fragment and recompress a little
		 * more to check whether it can be filled up.  Fix the fragment
		 * if succeeds.  Otherwise, just drop it and go on packing.
		 */
		if (may_packing && len == e->length &&
		    (compressedsize & (blksz - 1)) &&
		    ctx->tail < EROFS_COMPR_QUEUE_SZ) {
			fctx->pclustersize = roundup(compressedsize, blksz);
			goto fix_dedupedfrag;
		}

		if (may_inline && len == e->length)
			tryrecompress_trailing(ctx, h, ctx->queue + ctx->head,
					&e->length, dst, &compressedsize);

		e->compressedblks = BLK_ROUND_UP(sbi, compressedsize);
		DBG_BUGON(e->compressedblks * blksz >= e->length);

		padding = 0;
		tailused = compressedsize & (blksz - 1);
		if (tailused)
			padding = blksz - tailused;

		/* zero out garbage trailing data for non-0padding */
		if (!erofs_sb_has_lz4_0padding(sbi)) {
			memset(dst + compressedsize, 0, padding);
			padding = 0;
		}

		/* write compressed data */
		if (ctx->tmpfile) {
			erofs_dbg("Writing %u compressed data to tmpfile of %u blocks",
				  e->length, e->compressedblks);

			ret = fwrite(dst - padding, erofs_blksiz(sbi),
				     e->compressedblks, ctx->tmpfile);
			if (ret != e->compressedblks)
				return -EIO;
			fflush(ctx->tmpfile);
		} else {
			erofs_dbg("Writing %u compressed data to %u of %u blocks",
				  e->length, ctx->blkaddr, e->compressedblks);

			ret = blk_write(sbi, dst - padding, ctx->blkaddr,
					e->compressedblks);
			if (ret)
				return ret;
		}
		e->raw = false;
		may_inline = false;
		may_packing = false;
	}
	e->partial = false;
	e->blkaddr = ctx->blkaddr;
	if (!may_inline && !may_packing && !is_packed_inode)
		(void)z_erofs_dedupe_insert(e, ctx->queue + ctx->head);
	ctx->blkaddr += e->compressedblks;
	ctx->head += e->length;
	return 0;

fix_dedupedfrag:
	DBG_BUGON(!inode->fragment_size);
	ctx->remaining += inode->fragment_size;
	fctx->fix_dedupedfrag = true;
	return 1;
}

static int z_erofs_compress_one(struct z_erofs_vle_compress_ctx *ctx)
{
	unsigned int len = ctx->tail - ctx->head;
	struct z_erofs_extent_item *ei;
	struct z_erofs_file_compress_ctx *fctx = ctx->fctx;

	while (len) {
		int ret = z_erofs_compress_dedupe(ctx, &len);

		if (ret > 0)
			break;
		else if (ret < 0)
			return ret;

		DBG_BUGON(ctx->pivot);
		ei = malloc(sizeof(*ei));
		if (!ei)
			return -ENOMEM;

		init_list_head(&ei->list);
		ret = __z_erofs_compress_one(ctx, &ei->e);
		if (ret) {
			free(ei);
			if (ret > 0)
				break;		/* need more data */
			return ret;
		}

		len -= ei->e.length;
		ctx->pivot = ei;
		if (fctx->fix_dedupedfrag && !fctx->fragemitted &&
		    z_erofs_fixup_deduped_fragment(ctx, len))
			break;

		if (z_erofs_need_refill(ctx))
			break;
	}
	return 0;
}

struct z_erofs_compressindex_vec {
	union {
		erofs_blk_t blkaddr;
		u16 delta[2];
	} u;
	u16 clusterofs;
	u8  clustertype;
};

static void *parse_legacy_indexes(struct z_erofs_compressindex_vec *cv,
				  unsigned int nr, void *metacur)
{
	struct z_erofs_lcluster_index *const db = metacur;
	unsigned int i;

	for (i = 0; i < nr; ++i, ++cv) {
		struct z_erofs_lcluster_index *const di = db + i;
		const unsigned int advise = le16_to_cpu(di->di_advise);

		cv->clustertype = (advise >> Z_EROFS_LI_LCLUSTER_TYPE_BIT) &
			((1 << Z_EROFS_LI_LCLUSTER_TYPE_BITS) - 1);
		cv->clusterofs = le16_to_cpu(di->di_clusterofs);

		if (cv->clustertype == Z_EROFS_LCLUSTER_TYPE_NONHEAD) {
			cv->u.delta[0] = le16_to_cpu(di->di_u.delta[0]);
			cv->u.delta[1] = le16_to_cpu(di->di_u.delta[1]);
		} else {
			cv->u.blkaddr = le32_to_cpu(di->di_u.blkaddr);
		}
	}
	return db + nr;
}

static void *write_compacted_indexes(u8 *out,
				     struct z_erofs_compressindex_vec *cv,
				     erofs_blk_t *blkaddr_ret,
				     unsigned int destsize,
				     unsigned int lclusterbits,
				     bool final, bool *dummy_head,
				     bool update_blkaddr)
{
	unsigned int vcnt, lobits, encodebits, pos, i, cblks;
	erofs_blk_t blkaddr;

	if (destsize == 4)
		vcnt = 2;
	else if (destsize == 2 && lclusterbits <= 12)
		vcnt = 16;
	else
		return ERR_PTR(-EINVAL);
	lobits = max(lclusterbits, ilog2(Z_EROFS_LI_D0_CBLKCNT) + 1U);
	encodebits = (vcnt * destsize * 8 - 32) / vcnt;
	blkaddr = *blkaddr_ret;

	pos = 0;
	for (i = 0; i < vcnt; ++i) {
		unsigned int offset, v;
		u8 ch, rem;

		if (cv[i].clustertype == Z_EROFS_LCLUSTER_TYPE_NONHEAD) {
			if (cv[i].u.delta[0] & Z_EROFS_LI_D0_CBLKCNT) {
				cblks = cv[i].u.delta[0] & ~Z_EROFS_LI_D0_CBLKCNT;
				offset = cv[i].u.delta[0];
				blkaddr += cblks;
				*dummy_head = false;
			} else if (i + 1 == vcnt) {
				offset = min_t(u16, cv[i].u.delta[1],
						(1 << lobits) - 1);
			} else {
				offset = cv[i].u.delta[0];
			}
		} else {
			offset = cv[i].clusterofs;
			if (*dummy_head) {
				++blkaddr;
				if (update_blkaddr)
					*blkaddr_ret = blkaddr;
			}
			*dummy_head = true;
			update_blkaddr = false;

			if (cv[i].u.blkaddr != blkaddr) {
				if (i + 1 != vcnt)
					DBG_BUGON(!final);
				DBG_BUGON(cv[i].u.blkaddr);
			}
		}
		v = (cv[i].clustertype << lobits) | offset;
		rem = pos & 7;
		ch = out[pos / 8] & ((1 << rem) - 1);
		out[pos / 8] = (v << rem) | ch;
		out[pos / 8 + 1] = v >> (8 - rem);
		out[pos / 8 + 2] = v >> (16 - rem);
		pos += encodebits;
	}
	DBG_BUGON(destsize * vcnt * 8 != pos + 32);
	*(__le32 *)(out + destsize * vcnt - 4) = cpu_to_le32(*blkaddr_ret);
	*blkaddr_ret = blkaddr;
	return out + destsize * vcnt;
}

int z_erofs_convert_to_compacted_format(struct erofs_inode *inode,
					erofs_blk_t blkaddr,
					unsigned int legacymetasize,
					void *compressmeta)
{
	const unsigned int mpos = roundup(inode->inode_isize +
					  inode->xattr_isize, 8) +
				  sizeof(struct z_erofs_map_header);
	const unsigned int totalidx = (legacymetasize -
			Z_EROFS_LEGACY_MAP_HEADER_SIZE) /
				sizeof(struct z_erofs_lcluster_index);
	const unsigned int logical_clusterbits = inode->z_logical_clusterbits;
	u8 *out, *in;
	struct z_erofs_compressindex_vec cv[16];
	struct erofs_sb_info *sbi = inode->sbi;
	/* # of 8-byte units so that it can be aligned with 32 bytes */
	unsigned int compacted_4b_initial, compacted_4b_end;
	unsigned int compacted_2b;
	bool dummy_head;
	bool big_pcluster = erofs_sb_has_big_pcluster(sbi);

	if (logical_clusterbits < sbi->blkszbits)
		return -EINVAL;
	if (logical_clusterbits > 14) {
		erofs_err("compact format is unsupported for lcluster size %u",
			  1 << logical_clusterbits);
		return -EOPNOTSUPP;
	}

	if (inode->z_advise & Z_EROFS_ADVISE_COMPACTED_2B) {
		if (logical_clusterbits > 12) {
			erofs_err("compact 2B is unsupported for lcluster size %u",
				  1 << logical_clusterbits);
			return -EINVAL;
		}

		compacted_4b_initial = (32 - mpos % 32) / 4;
		if (compacted_4b_initial == 32 / 4)
			compacted_4b_initial = 0;

		if (compacted_4b_initial > totalidx) {
			compacted_4b_initial = compacted_2b = 0;
			compacted_4b_end = totalidx;
		} else {
			compacted_2b = rounddown(totalidx -
						 compacted_4b_initial, 16);
			compacted_4b_end = totalidx - compacted_4b_initial -
					   compacted_2b;
		}
	} else {
		compacted_2b = compacted_4b_initial = 0;
		compacted_4b_end = totalidx;
	}

	out = in = compressmeta;

	out += sizeof(struct z_erofs_map_header);
	in += Z_EROFS_LEGACY_MAP_HEADER_SIZE;

	dummy_head = false;
	/* prior to bigpcluster, blkaddr was bumped up once coming into HEAD */
	if (!big_pcluster) {
		--blkaddr;
		dummy_head = true;
	}

	/* generate compacted_4b_initial */
	while (compacted_4b_initial) {
		in = parse_legacy_indexes(cv, 2, in);
		out = write_compacted_indexes(out, cv, &blkaddr,
					      4, logical_clusterbits, false,
					      &dummy_head, big_pcluster);
		compacted_4b_initial -= 2;
	}
	DBG_BUGON(compacted_4b_initial);

	/* generate compacted_2b */
	while (compacted_2b) {
		in = parse_legacy_indexes(cv, 16, in);
		out = write_compacted_indexes(out, cv, &blkaddr,
					      2, logical_clusterbits, false,
					      &dummy_head, big_pcluster);
		compacted_2b -= 16;
	}
	DBG_BUGON(compacted_2b);

	/* generate compacted_4b_end */
	while (compacted_4b_end > 1) {
		in = parse_legacy_indexes(cv, 2, in);
		out = write_compacted_indexes(out, cv, &blkaddr,
					      4, logical_clusterbits, false,
					      &dummy_head, big_pcluster);
		compacted_4b_end -= 2;
	}

	/* generate final compacted_4b_end if needed */
	if (compacted_4b_end) {
		memset(cv, 0, sizeof(cv));
		in = parse_legacy_indexes(cv, 1, in);
		out = write_compacted_indexes(out, cv, &blkaddr,
					      4, logical_clusterbits, true,
					      &dummy_head, big_pcluster);
	}
	inode->extent_isize = out - (u8 *)compressmeta;
	return 0;
}

static void z_erofs_write_mapheader(struct erofs_inode *inode,
				    void *compressmeta)
{
	struct erofs_sb_info *sbi = inode->sbi;
	struct z_erofs_map_header h = {
		.h_advise = cpu_to_le16(inode->z_advise),
		.h_algorithmtype = inode->z_algorithmtype[1] << 4 |
				   inode->z_algorithmtype[0],
		/* lclustersize */
		.h_clusterbits = inode->z_logical_clusterbits - sbi->blkszbits,
	};

	if (inode->z_advise & Z_EROFS_ADVISE_FRAGMENT_PCLUSTER)
		h.h_fragmentoff = cpu_to_le32(inode->fragmentoff);
	else
		h.h_idata_size = cpu_to_le16(inode->idata_size);

	memset(compressmeta, 0, Z_EROFS_LEGACY_MAP_HEADER_SIZE);
	/* write out map header */
	memcpy(compressmeta, &h, sizeof(struct z_erofs_map_header));
}

void z_erofs_drop_inline_pcluster(struct erofs_inode *inode)
{
	struct erofs_sb_info *sbi = inode->sbi;
	const unsigned int type = Z_EROFS_LCLUSTER_TYPE_PLAIN;
	struct z_erofs_map_header *h = inode->compressmeta;

	h->h_advise = cpu_to_le16(le16_to_cpu(h->h_advise) &
				  ~Z_EROFS_ADVISE_INLINE_PCLUSTER);
	h->h_idata_size = 0;
	if (!inode->eof_tailraw)
		return;
	DBG_BUGON(inode->compressed_idata != true);

	/* patch the EOF lcluster to uncompressed type first */
	if (inode->datalayout == EROFS_INODE_COMPRESSED_FULL) {
		struct z_erofs_lcluster_index *di =
			(inode->compressmeta + inode->extent_isize) -
			sizeof(struct z_erofs_lcluster_index);
		__le16 advise =
			cpu_to_le16(type << Z_EROFS_LI_LCLUSTER_TYPE_BIT);

		di->di_advise = advise;
	} else if (inode->datalayout == EROFS_INODE_COMPRESSED_COMPACT) {
		/* handle the last compacted 4B pack */
		unsigned int eofs, base, pos, v, lo;
		u8 *out;

		eofs = inode->extent_isize -
			(4 << (BLK_ROUND_UP(sbi, inode->i_size) & 1));
		base = round_down(eofs, 8);
		pos = 16 /* encodebits */ * ((eofs - base) / 4);
		out = inode->compressmeta + base;
		lo = erofs_blkoff(sbi, get_unaligned_le32(out + pos / 8));
		v = (type << sbi->blkszbits) | lo;
		out[pos / 8] = v & 0xff;
		out[pos / 8 + 1] = v >> 8;
	} else {
		DBG_BUGON(1);
		return;
	}
	free(inode->idata);
	/* replace idata with prepared uncompressed data */
	inode->idata = inode->eof_tailraw;
	inode->idata_size = inode->eof_tailrawsize;
	inode->compressed_idata = false;
	inode->eof_tailraw = NULL;
}

int z_erofs_compress_file(struct z_erofs_vle_compress_ctx *ctx, u64 offset,
			  erofs_blk_t blkaddr)
{
	struct z_erofs_file_compress_ctx *fctx = ctx->fctx;
	struct erofs_inode *inode = fctx->inode;
	int ret = 0;

	while (ctx->remaining) {
		const u64 rx = min_t(u64, ctx->remaining,
				     EROFS_COMPR_QUEUE_SZ - ctx->tail);

		ret = pread(fctx->fd, ctx->queue + ctx->tail, rx, offset);
		if (ret != rx)
			return -errno;
		ctx->remaining -= rx;
		ctx->tail += rx;
		offset += rx;

		ret = z_erofs_compress_one(ctx);
		if (ret)
			return ret;
	}
	DBG_BUGON(ctx->head != ctx->tail);

	ctx->compressed_blocks = ctx->blkaddr - blkaddr;
	DBG_BUGON(ctx->compressed_blocks < !!inode->idata_size);
	ctx->compressed_blocks -= !!inode->idata_size;

	if (ctx->pivot) {
		z_erofs_commit_extent(ctx, ctx->pivot);
		ctx->pivot = NULL;
	}

	return 0;
}

#ifdef EROFS_MT_ENABLED
int z_erofs_mt_private_init(struct erofs_sb_info *sbi,
			    struct erofs_compress_wq_private *priv,
			    unsigned int alg_id, char *alg_name,
			    unsigned int comp_level, unsigned int dict_size)
{
	struct erofs_compress_cfg *lc;
	int ret;

	if (unlikely(!priv->init)) {
		priv->init = true;

		priv->queue = malloc(EROFS_COMPR_QUEUE_SZ);
		if (!priv->queue)
			return -ENOMEM;

		priv->destbuf = calloc(1, EROFS_CONFIG_COMPR_MAX_SZ +
						  EROFS_MAX_BLOCK_SIZE);
		if (!priv->destbuf)
			return -ENOMEM;

		priv->ccfg = calloc(EROFS_MAX_COMPR_CFGS,
				    sizeof(struct erofs_compress_cfg));
		if (!priv->ccfg)
			return -ENOMEM;
#ifdef USE_PER_WORKER_TMPFILE
#ifndef HAVE_TMPFILE64
		priv->tmpfile = tmpfile();
#else
		priv->tmpfile = tmpfile64();
#endif
		if (!priv->tmpfile)
			return -errno;
#endif
	}

	lc = &priv->ccfg[alg_id];
	if (!lc->enable) {
		lc->enable = true;
		lc->algorithmtype = alg_id;

		ret = erofs_compressor_init(sbi, &lc->handle, alg_name,
					    comp_level, dict_size);
		if (ret)
			return ret;
	}

	return 0;
}

void z_erofs_mt_private_fini(void *private)
{
	struct erofs_compress_wq_private *priv = private;
	int i;

	if (priv->init) {
		for (i = 0; i < EROFS_MAX_COMPR_CFGS; i++) {
			if (priv->ccfg[i].enable)
				erofs_compressor_exit(&priv->ccfg[i].handle);
		}
		free(priv->ccfg);
		free(priv->destbuf);
		free(priv->queue);
#ifdef USE_PER_WORKER_TMPFILE
		fclose(priv->tmpfile);
#endif
		priv->init = false;
	}
}

void z_erofs_mt_work(struct erofs_work *work)
{
	struct erofs_compress_work *cwork = (struct erofs_compress_work *)work;
	struct z_erofs_vle_compress_ctx *ctx = &cwork->ctx;
	struct erofs_compress_wq_private *priv = work->priv;
	erofs_blk_t blkaddr = ctx->blkaddr;
	u64 offset = ctx->seg_idx * cfg.c_segment_size;
	int ret = 0;

	ret = z_erofs_mt_private_init(ctx->fctx->inode->sbi, priv,
				      cwork->alg_id, cwork->alg_name,
				      cwork->comp_level, cwork->dict_size);
	if (ret)
		goto out;

	ctx->queue = priv->queue;
	ctx->destbuf = priv->destbuf;
	ctx->chandle = &priv->ccfg[cwork->alg_id].handle;
#ifdef USE_PER_WORKER_TMPFILE
	ctx->tmpfile = priv->tmpfile;
	ctx->tmpfile_off = ftell(ctx->tmpfile);
	if (ctx->tmpfile_off == -1) {
		ret = -errno;
		goto out;
	}
#else
#ifdef HAVE_TMPFILE64
	ctx->tmpfile = tmpfile64();
#else
	ctx->tmpfile = tmpfile();
#endif
	if (!ctx->tmpfile) {
		ret = -errno;
		goto out;
	}
	ctx->tmpfile_off = 0;
#endif

	ret = z_erofs_compress_file(ctx, offset, blkaddr);
	if (ret)
		goto out;

	fflush(ctx->tmpfile);

out:
	cwork->ret = ret;
	pthread_mutex_lock(&z_erofs_mt_ctrl.mutex);
	++z_erofs_mt_ctrl.nfini;
	pthread_cond_signal(&z_erofs_mt_ctrl.cond);
	pthread_mutex_unlock(&z_erofs_mt_ctrl.mutex);
}

int z_erofs_mt_merge(struct erofs_compress_work *cur, erofs_blk_t blkaddr,
		     erofs_blk_t *compressed_blocks)
{
	struct z_erofs_vle_compress_ctx *ctx, *listhead = NULL;
	struct erofs_sb_info *sbi = cur->ctx.fctx->inode->sbi;
	struct erofs_compress_work *tmp;
	char *memblock = NULL;
	size_t size = 0;
	int ret = 0, lret;

	while (cur != NULL) {
		ctx = &cur->ctx;

		if (!listhead)
			listhead = ctx;
		else
			list_splice_tail(&ctx->extents, &listhead->extents);

		if (cur->ret != 0) {
			if (!ret)
				ret = cur->ret;
			goto out;
		}

		size = ctx->compressed_blocks * erofs_blksiz(sbi);
		memblock = realloc(memblock, size);
		if (!memblock) {
			if (!ret)
				ret = -ENOMEM;
			goto out;
		}

		lret = pread(fileno(ctx->tmpfile), memblock, size,
			     ctx->tmpfile_off);
		if (lret != size) {
			if (!ret)
				ret = errno ? -errno : -EIO;
			goto out;
		}

#ifdef USE_PER_WORKER_TMPFILE
		lret = fallocate(fileno(ctx->tmpfile),
				 FALLOC_FL_PUNCH_HOLE | FALLOC_FL_KEEP_SIZE,
				 ctx->tmpfile_off, size);
		if (lret) {
			if (!ret)
				ret = -errno;
			goto out;
		}
#endif

		lret = blk_write(sbi, memblock, blkaddr + *compressed_blocks,
				 ctx->compressed_blocks);
		if (lret) {
			if (!ret)
				ret = lret;
			goto out;
		}
		*compressed_blocks += ctx->compressed_blocks;

out:
#ifndef USE_PER_WORKER_TMPFILE
		fclose(ctx->tmpfile);
#endif
		tmp = cur->next;
		cur->next = z_erofs_mt_ctrl.idle;
		z_erofs_mt_ctrl.idle = cur;
		cur = tmp;
	}

	free(memblock);

	return ret;
}

int z_erofs_mt_compress(struct z_erofs_vle_compress_ctx *ctx,
			struct z_erofs_write_index_ctx *ictx,
			struct erofs_compress_cfg *ccfg,
			erofs_blk_t *compressed_blocks)
{
	struct erofs_compress_work *work, *head = NULL, **last = &head;
	struct z_erofs_file_compress_ctx *fctx = ctx->fctx;
	struct erofs_inode *inode = fctx->inode;
	erofs_blk_t blkaddr = ctx->blkaddr;
	int nsegs = DIV_ROUND_UP(inode->i_size, cfg.c_segment_size);

	z_erofs_mt_ctrl.nfini = 0;

	for (int i = 0; i < nsegs; i++) {
		if (z_erofs_mt_ctrl.idle) {
			work = z_erofs_mt_ctrl.idle;
			z_erofs_mt_ctrl.idle = work->next;
			work->next = NULL;
		} else {
			work = calloc(1, sizeof(*work));
			if (!work)
				return -ENOMEM;
		}
		*last = work;
		last = &work->next;

		memset(&work->ctx, 0, sizeof(work->ctx));
		if (i == nsegs - 1)
			work->ctx.remaining = inode->i_size -
					      inode->fragment_size -
					      i * cfg.c_segment_size;
		else
			work->ctx.remaining = cfg.c_segment_size;
		work->ctx.seg_num = nsegs;
		work->ctx.seg_idx = i;
		work->ctx.blkaddr = blkaddr;
		init_list_head(&work->ctx.extents);
		work->ctx.fctx = fctx;

		work->alg_id = ccfg->handle.alg->id;
		work->alg_name = ccfg->handle.alg->name;
		work->comp_level = ccfg->handle.compression_level;
		work->dict_size = ccfg->handle.dict_size;

		work->work.func = z_erofs_mt_work;

		erofs_queue_work(&z_erofs_mt_ctrl.wq, &work->work);
	}

	pthread_mutex_lock(&z_erofs_mt_ctrl.mutex);
	while (z_erofs_mt_ctrl.nfini != nsegs)
		pthread_cond_wait(&z_erofs_mt_ctrl.cond,
				  &z_erofs_mt_ctrl.mutex);
	pthread_mutex_unlock(&z_erofs_mt_ctrl.mutex);

	ictx->extents = &head->ctx.extents;

	return z_erofs_mt_merge(head, blkaddr, compressed_blocks);
}
#endif

int erofs_write_compressed_file(struct erofs_inode *inode, int fd)
{
	struct erofs_buffer_head *bh;
	static struct z_erofs_file_compress_ctx fctx;
	static struct z_erofs_vle_compress_ctx ctx;
	static struct z_erofs_write_index_ctx ictx;
	struct erofs_compress_cfg *ccfg;
	erofs_blk_t blkaddr, compressed_blocks = 0;
	unsigned int legacymetasize;
	int ret;
	struct erofs_sb_info *sbi = inode->sbi;
	u8 *compressmeta = malloc(BLK_ROUND_UP(sbi, inode->i_size) *
				  sizeof(struct z_erofs_lcluster_index) +
				  Z_EROFS_LEGACY_MAP_HEADER_SIZE);

	if (!compressmeta)
		return -ENOMEM;

	/* allocate main data buffer */
	bh = erofs_balloc(DATA, 0, 0, 0);
	if (IS_ERR(bh)) {
		ret = PTR_ERR(bh);
		goto err_free_meta;
	}

	/* initialize per-file compression setting */
	inode->z_advise = 0;
	inode->z_logical_clusterbits = sbi->blkszbits;
	if (!cfg.c_legacy_compress && inode->z_logical_clusterbits <= 14) {
		if (inode->z_logical_clusterbits <= 12)
			inode->z_advise |= Z_EROFS_ADVISE_COMPACTED_2B;
		inode->datalayout = EROFS_INODE_COMPRESSED_COMPACT;
	} else {
		inode->datalayout = EROFS_INODE_COMPRESSED_FULL;
	}

	if (erofs_sb_has_big_pcluster(sbi)) {
		inode->z_advise |= Z_EROFS_ADVISE_BIG_PCLUSTER_1;
		if (inode->datalayout == EROFS_INODE_COMPRESSED_COMPACT)
			inode->z_advise |= Z_EROFS_ADVISE_BIG_PCLUSTER_2;
	}
	if (cfg.c_fragments && !cfg.c_dedupe)
		inode->z_advise |= Z_EROFS_ADVISE_INTERLACED_PCLUSTER;

#ifndef NDEBUG
	if (cfg.c_random_algorithms) {
		while (1) {
			inode->z_algorithmtype[0] =
				rand() % EROFS_MAX_COMPR_CFGS;
			if (erofs_ccfg[inode->z_algorithmtype[0]].enable)
				break;
		}
	}
#endif
	ccfg = &erofs_ccfg[inode->z_algorithmtype[0]];
	inode->z_algorithmtype[0] = ccfg[0].algorithmtype;
	inode->z_algorithmtype[1] = 0;

	inode->idata_size = 0;
	inode->fragment_size = 0;

	/*
	 * Handle tails in advance to avoid writing duplicated
	 * parts into the packed inode.
	 */
	if (cfg.c_fragments && !erofs_is_packed_inode(inode)) {
		ret = z_erofs_fragments_dedupe(inode, fd, &fctx.tof_chksum);
		if (ret < 0)
			goto err_bdrop;
	}

	blkaddr = erofs_mapbh(bh->block);	/* start_blkaddr */

	fctx.inode = inode;
	fctx.fd = fd;
	fctx.fix_dedupedfrag = false;
	fctx.fragemitted = false;
	fctx.pclustersize = z_erofs_get_max_pclustersize(inode);

	memset(&ctx, 0, sizeof(ctx));
	ctx.fctx = &fctx;
	ctx.blkaddr = blkaddr;
	init_list_head(&ctx.extents);

	if (z_erofs_mt_enabled) {
#ifdef EROFS_MT_ENABLED
		if (inode->i_size <= cfg.c_segment_size)
			goto single_thread_comp;

		ret = z_erofs_mt_compress(&ctx, &ictx, ccfg,
					  &compressed_blocks);
		if (ret)
			goto err_free_idata;
#endif
	} else {
		if (cfg.c_all_fragments && !erofs_is_packed_inode(inode) &&
		    !inode->fragment_size) {
			ret = z_erofs_pack_file_from_fd(inode, fd,
							fctx.tof_chksum);
			if (ret)
				goto err_free_idata;

			ictx.extents = &ctx.extents;
		} else {
#ifdef EROFS_MT_ENABLED
single_thread_comp:
#endif
			ctx.queue = z_erofs_global_queue;
			ctx.destbuf = NULL;
			ctx.chandle = &ccfg->handle;
			ctx.remaining = inode->i_size - inode->fragment_size;
			ctx.seg_num = 1;
			ctx.seg_idx = 0;

			ret = z_erofs_compress_file(&ctx, 0, blkaddr);
			if (ret)
				goto err_free_idata;

			compressed_blocks = ctx.compressed_blocks;
			ictx.extents = &ctx.extents;
		}

		/* generate an extent for the deduplicated fragment */
		if (inode->fragment_size && !fctx.fragemitted) {
			struct z_erofs_extent_item *ei;

			ei = malloc(sizeof(*ei));
			if (!ei) {
				ret = -ENOMEM;
				goto err_free_idata;
			}

			ei->e = (struct z_erofs_inmem_extent){
				.length = inode->fragment_size,
				.compressedblks = 0,
				.raw = false,
				.partial = false,
				.blkaddr = ctx.blkaddr,
			};
			init_list_head(&ei->list);
			z_erofs_commit_extent(&ctx, ei);
		}
		z_erofs_fragments_commit(inode);
	}

	ictx.inode = inode;
	ictx.blkaddr = blkaddr;
	ictx.blkoff = 0;
	ictx.clusterofs = 0;
	ictx.metacur = compressmeta + Z_EROFS_LEGACY_MAP_HEADER_SIZE;

	z_erofs_write_indexes(&ictx);
	legacymetasize = ictx.metacur - compressmeta;
	/* estimate if data compression saves space or not */
	if (!inode->fragment_size &&
	    compressed_blocks * erofs_blksiz(sbi) + inode->idata_size +
	    legacymetasize >= inode->i_size) {
		z_erofs_dedupe_commit(true);
		ret = -ENOSPC;
		goto err_free_idata;
	}
	z_erofs_dedupe_commit(false);
	z_erofs_write_mapheader(inode, compressmeta);

	if (!fctx.fragemitted)
		sbi->saved_by_deduplication += inode->fragment_size;

	/* if the entire file is a fragment, a simplified form is used. */
	if (inode->i_size <= inode->fragment_size) {
		DBG_BUGON(inode->i_size < inode->fragment_size);
		DBG_BUGON(inode->fragmentoff >> 63);
		*(__le64 *)compressmeta =
			cpu_to_le64(inode->fragmentoff | 1ULL << 63);
		inode->datalayout = EROFS_INODE_COMPRESSED_FULL;
		legacymetasize = Z_EROFS_LEGACY_MAP_HEADER_SIZE;
	}

	if (compressed_blocks) {
		ret = erofs_bh_balloon(bh, erofs_pos(sbi, compressed_blocks));
		DBG_BUGON(ret != erofs_blksiz(sbi));
	} else {
		if (!cfg.c_fragments && !cfg.c_dedupe)
			DBG_BUGON(!inode->idata_size);
	}

	erofs_info("compressed %s (%llu bytes) into %u blocks",
		   inode->i_srcpath, (unsigned long long)inode->i_size,
		   compressed_blocks);

	if (inode->idata_size) {
		bh->op = &erofs_skip_write_bhops;
		inode->bh_data = bh;
	} else {
		erofs_bdrop(bh, false);
	}

	inode->u.i_blocks = compressed_blocks;

	if (inode->datalayout == EROFS_INODE_COMPRESSED_FULL) {
		inode->extent_isize = legacymetasize;
	} else {
		ret = z_erofs_convert_to_compacted_format(inode, blkaddr,
							  legacymetasize,
							  compressmeta);
		DBG_BUGON(ret);
	}
	inode->compressmeta = compressmeta;
	if (!erofs_is_packed_inode(inode))
		erofs_droid_blocklist_write(inode, blkaddr, compressed_blocks);
	return 0;

err_free_idata:
	if (inode->idata) {
		free(inode->idata);
		inode->idata = NULL;
	}
err_bdrop:
	erofs_bdrop(bh, true);	/* revoke buffer */
err_free_meta:
	free(compressmeta);
	inode->compressmeta = NULL;
	return ret;
}

static int z_erofs_build_compr_cfgs(struct erofs_sb_info *sbi,
				    struct erofs_buffer_head *sb_bh,
				    u32 *max_dict_size)
{
	struct erofs_buffer_head *bh = sb_bh;
	int ret = 0;

	if (sbi->available_compr_algs & (1 << Z_EROFS_COMPRESSION_LZ4)) {
		struct {
			__le16 size;
			struct z_erofs_lz4_cfgs lz4;
		} __packed lz4alg = {
			.size = cpu_to_le16(sizeof(struct z_erofs_lz4_cfgs)),
			.lz4 = {
				.max_distance =
					cpu_to_le16(sbi->lz4_max_distance),
				.max_pclusterblks = cfg.c_pclusterblks_max,
			}
		};

		bh = erofs_battach(bh, META, sizeof(lz4alg));
		if (IS_ERR(bh)) {
			DBG_BUGON(1);
			return PTR_ERR(bh);
		}
		erofs_mapbh(bh->block);
		ret = dev_write(sbi, &lz4alg, erofs_btell(bh, false),
				sizeof(lz4alg));
		bh->op = &erofs_drop_directly_bhops;
	}
#ifdef HAVE_LIBLZMA
	if (sbi->available_compr_algs & (1 << Z_EROFS_COMPRESSION_LZMA)) {
		struct {
			__le16 size;
			struct z_erofs_lzma_cfgs lzma;
		} __packed lzmaalg = {
			.size = cpu_to_le16(sizeof(struct z_erofs_lzma_cfgs)),
			.lzma = {
				.dict_size = cpu_to_le32(
					max_dict_size
						[Z_EROFS_COMPRESSION_LZMA]),
			}
		};

		bh = erofs_battach(bh, META, sizeof(lzmaalg));
		if (IS_ERR(bh)) {
			DBG_BUGON(1);
			return PTR_ERR(bh);
		}
		erofs_mapbh(bh->block);
		ret = dev_write(sbi, &lzmaalg, erofs_btell(bh, false),
				sizeof(lzmaalg));
		bh->op = &erofs_drop_directly_bhops;
	}
#endif
	if (sbi->available_compr_algs & (1 << Z_EROFS_COMPRESSION_DEFLATE)) {
		struct {
			__le16 size;
			struct z_erofs_deflate_cfgs z;
		} __packed zalg = {
			.size = cpu_to_le16(sizeof(struct z_erofs_deflate_cfgs)),
			.z = {
				.windowbits = cpu_to_le32(ilog2(
					max_dict_size
						[Z_EROFS_COMPRESSION_DEFLATE])),
			}
		};

		bh = erofs_battach(bh, META, sizeof(zalg));
		if (IS_ERR(bh)) {
			DBG_BUGON(1);
			return PTR_ERR(bh);
		}
		erofs_mapbh(bh->block);
		ret = dev_write(sbi, &zalg, erofs_btell(bh, false),
				sizeof(zalg));
		bh->op = &erofs_drop_directly_bhops;
	}
	return ret;
}

int z_erofs_compress_init(struct erofs_sb_info *sbi, struct erofs_buffer_head *sb_bh)
{
	int i, ret, id;
	u32 max_dict_size[Z_EROFS_COMPRESSION_MAX] = {};

	for (i = 0; cfg.c_compr_opts[i].alg; ++i) {
		struct erofs_compress *c = &erofs_ccfg[i].handle;

		ret = erofs_compressor_init(sbi, c, cfg.c_compr_opts[i].alg,
					    cfg.c_compr_opts[i].level,
					    cfg.c_compr_opts[i].dict_size);
		if (ret)
			return ret;

		id = z_erofs_get_compress_algorithm_id(c);
		erofs_ccfg[i].algorithmtype = id;
		erofs_ccfg[i].enable = true;
		sbi->available_compr_algs |= 1 << erofs_ccfg[i].algorithmtype;
		if (erofs_ccfg[i].algorithmtype != Z_EROFS_COMPRESSION_LZ4)
			erofs_sb_set_compr_cfgs(sbi);
		if (c->dict_size > max_dict_size[id])
			max_dict_size[id] = c->dict_size;
	}

	/*
	 * if primary algorithm is empty (e.g. compression off),
	 * clear 0PADDING feature for old kernel compatibility.
	 */
	if (!cfg.c_compr_opts[0].alg ||
	    (cfg.c_legacy_compress &&
	     !strncmp(cfg.c_compr_opts[0].alg, "lz4", 3)))
		erofs_sb_clear_lz4_0padding(sbi);

	if (!cfg.c_compr_opts[0].alg)
		return 0;

	/*
	 * if big pcluster is enabled, an extra CBLKCNT lcluster index needs
	 * to be loaded in order to get those compressed block counts.
	 */
	if (cfg.c_pclusterblks_max > 1) {
		if (cfg.c_pclusterblks_max >
		    Z_EROFS_PCLUSTER_MAX_SIZE / erofs_blksiz(sbi)) {
			erofs_err("unsupported clusterblks %u (too large)",
				  cfg.c_pclusterblks_max);
			return -EINVAL;
		}
		erofs_sb_set_big_pcluster(sbi);
	}
	if (cfg.c_pclusterblks_packed > cfg.c_pclusterblks_max) {
		erofs_err("invalid physical cluster size for the packed file");
		return -EINVAL;
	}

	if (erofs_sb_has_compr_cfgs(sbi)) {
		ret = z_erofs_build_compr_cfgs(sbi, sb_bh, max_dict_size);
		if (ret)
			return ret;
	}

#ifdef EROFS_MT_ENABLED
	if (cfg.c_mt_workers == 1) {
		z_erofs_mt_enabled = false;
	} else {
		pthread_mutex_init(&z_erofs_mt_ctrl.mutex, NULL);
		pthread_cond_init(&z_erofs_mt_ctrl.cond, NULL);
		ret = erofs_alloc_workqueue(
			&z_erofs_mt_ctrl.wq, cfg.c_mt_workers,
			cfg.c_mt_workers << 2,
			sizeof(struct erofs_compress_wq_private),
			z_erofs_mt_private_fini);
		z_erofs_mt_enabled = !ret;
	}
#else
	mt_enabled = false;
#endif
	z_erofs_global_queue = malloc(EROFS_COMPR_QUEUE_SZ);
	if (!z_erofs_global_queue)
		return -ENOMEM;

	return 0;
}

int z_erofs_compress_exit(void)
{
	int i, ret;

	for (i = 0; cfg.c_compr_opts[i].alg; ++i) {
		ret = erofs_compressor_exit(&erofs_ccfg[i].handle);
		if (ret)
			return ret;
	}

	if (z_erofs_mt_enabled) {
#ifdef EROFS_MT_ENABLED
		ret = erofs_destroy_workqueue(&z_erofs_mt_ctrl.wq);
		if (ret)
			return ret;
		while (z_erofs_mt_ctrl.idle) {
			struct erofs_compress_work *tmp =
				z_erofs_mt_ctrl.idle->next;
			free(z_erofs_mt_ctrl.idle);
			z_erofs_mt_ctrl.idle = tmp;
		}
#endif
	}

	free(z_erofs_global_queue);

	return 0;
}
