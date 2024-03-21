// SPDX-License-Identifier: GPL-2.0+ OR Apache-2.0
/*
 * Copyright (C) 2018-2019 HUAWEI, Inc.
 *             http://www.huawei.com/
 * Created by Miao Xie <miaoxie@huawei.com>
 * with heavy changes by Gao Xiang <xiang@kernel.org>
 */
#ifndef _LARGEFILE64_SOURCE
#define _LARGEFILE64_SOURCE
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

struct z_erofs_compress_ictx {		/* inode context */
	struct erofs_inode *inode;
	unsigned int pclustersize;
	int fd;
	u64 fpos;

	u32 tof_chksum;
	bool fix_dedupedfrag;
	bool fragemitted;

	/* fields for write indexes */
	u8 *metacur;
	struct list_head extents;
	u16 clusterofs;
};

struct z_erofs_compress_sctx {		/* segment context */
	struct z_erofs_compress_ictx *ictx;

	u8 *queue;
	struct list_head extents;
	struct z_erofs_extent_item *pivot;

	struct erofs_compress *chandle;
	char *destbuf;

	unsigned int head, tail;
	erofs_off_t remaining;
	erofs_blk_t blkaddr;		/* pointing to the next blkaddr */
	u16 clusterofs;

	int seg_num, seg_idx;

	void *membuf;
	erofs_off_t memoff;
};

#ifdef EROFS_MT_ENABLED
struct erofs_compress_wq_tls {
	u8 *queue;
	char *destbuf;
	struct erofs_compress_cfg *ccfg;
};

struct erofs_compress_work {
	/* Note: struct erofs_work must be the first member */
	struct erofs_work work;
	struct z_erofs_compress_sctx ctx;
	struct erofs_compress_work *next;
	struct z_erofs_mt_file *mtfile_desc;

	unsigned int alg_id;
	char *alg_name;
	unsigned int comp_level;
	unsigned int dict_size;

	int errcode;
};

static struct {
	struct erofs_workqueue wq;
	struct erofs_compress_work *work_idle;
	pthread_mutex_t work_mutex;
	struct z_erofs_mt_file *file_idle;
	pthread_mutex_t file_mutex;
} z_erofs_mt_ctrl;
#endif

bool z_erofs_mt_enabled;

#define Z_EROFS_LEGACY_MAP_HEADER_SIZE	Z_EROFS_FULL_INDEX_ALIGN(0)

static void z_erofs_write_indexes_final(struct z_erofs_compress_ictx *ctx)
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

static void z_erofs_write_extent(struct z_erofs_compress_ictx *ctx,
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
		    !e->compressedblks)
			di.di_u.blkaddr = cpu_to_le32(inode->fragmentoff >> 32);
		else
			di.di_u.blkaddr = cpu_to_le32(e->blkaddr);
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
			    !e->compressedblks)
				di.di_u.blkaddr = cpu_to_le32(inode->fragmentoff >> 32);
			else
				di.di_u.blkaddr = cpu_to_le32(e->blkaddr);

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

static void z_erofs_write_indexes(struct z_erofs_compress_ictx *ctx)
{
	struct z_erofs_extent_item *ei, *n;

	ctx->clusterofs = 0;
	list_for_each_entry_safe(ei, n, &ctx->extents, list) {
		z_erofs_write_extent(ctx, &ei->e);

		list_del(&ei->list);
		free(ei);
	}
	z_erofs_write_indexes_final(ctx);
}

static bool z_erofs_need_refill(struct z_erofs_compress_sctx *ctx)
{
	const bool final = !ctx->remaining;
	unsigned int qh_aligned, qh_after;
	struct erofs_inode *inode = ctx->ictx->inode;

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

static void z_erofs_commit_extent(struct z_erofs_compress_sctx *ctx,
				  struct z_erofs_extent_item *ei)
{
	if (ei == &dummy_pivot)
		return;

	list_add_tail(&ei->list, &ctx->extents);
	ctx->clusterofs = (ctx->clusterofs + ei->e.length) &
			  (erofs_blksiz(ctx->ictx->inode->sbi) - 1);
}

static int z_erofs_compress_dedupe(struct z_erofs_compress_sctx *ctx,
				   unsigned int *len)
{
	struct erofs_inode *inode = ctx->ictx->inode;
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

static int write_uncompressed_extent(struct z_erofs_compress_sctx *ctx,
				     unsigned int len, char *dst)
{
	struct erofs_inode *inode = ctx->ictx->inode;
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

	if (ctx->membuf) {
		erofs_dbg("Writing %u uncompressed data of %s", count,
			  inode->i_srcpath);
		memcpy(ctx->membuf + ctx->memoff, dst, erofs_blksiz(sbi));
		ctx->memoff += erofs_blksiz(sbi);
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

static void tryrecompress_trailing(struct z_erofs_compress_sctx *ctx,
				   struct erofs_compress *ec,
				   void *in, unsigned int *insize,
				   void *out, unsigned int *compressedsize)
{
	struct erofs_sb_info *sbi = ctx->ictx->inode->sbi;
	static char tmp[Z_EROFS_PCLUSTER_MAX_SIZE];
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

static bool z_erofs_fixup_deduped_fragment(struct z_erofs_compress_sctx *ctx,
					   unsigned int len)
{
	struct z_erofs_compress_ictx *ictx = ctx->ictx;
	struct erofs_inode *inode = ictx->inode;
	struct erofs_sb_info *sbi = inode->sbi;
	const unsigned int newsize = ctx->remaining + len;

	DBG_BUGON(!inode->fragment_size);

	/* try to fix again if it gets larger (should be rare) */
	if (inode->fragment_size < newsize) {
		ictx->pclustersize = min_t(erofs_off_t,
				z_erofs_get_max_pclustersize(inode),
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

static int __z_erofs_compress_one(struct z_erofs_compress_sctx *ctx,
				  struct z_erofs_inmem_extent *e)
{
	static char g_dstbuf[EROFS_CONFIG_COMPR_MAX_SZ + EROFS_MAX_BLOCK_SIZE];
	char *dstbuf = ctx->destbuf ?: g_dstbuf;
	struct z_erofs_compress_ictx *ictx = ctx->ictx;
	struct erofs_inode *inode = ictx->inode;
	struct erofs_sb_info *sbi = inode->sbi;
	unsigned int blksz = erofs_blksiz(sbi);
	char *const dst = dstbuf + blksz;
	struct erofs_compress *const h = ctx->chandle;
	unsigned int len = ctx->tail - ctx->head;
	bool is_packed_inode = erofs_is_packed_inode(inode);
	bool tsg = (ctx->seg_idx + 1 >= ctx->seg_num), final = !ctx->remaining;
	bool may_packing = (cfg.c_fragments && tsg && final &&
			    !is_packed_inode && !z_erofs_mt_enabled);
	bool may_inline = (cfg.c_ztailpacking && tsg && final && !may_packing);
	unsigned int compressedsize;
	int ret;

	if (len <= ictx->pclustersize) {
		if (!final || !len)
			return 1;
		if (may_packing) {
			if (inode->fragment_size && !ictx->fix_dedupedfrag) {
				ictx->pclustersize = roundup(len, blksz);
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
				      &e->length, dst, ictx->pclustersize);
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
		   compressedsize < ictx->pclustersize &&
		   (!inode->fragment_size || ictx->fix_dedupedfrag)) {
frag_packing:
		ret = z_erofs_pack_fragments(inode, ctx->queue + ctx->head,
					     len, ictx->tof_chksum);
		if (ret < 0)
			return ret;
		e->compressedblks = 0; /* indicate a fragment */
		e->raw = false;
		ictx->fragemitted = true;
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
		    ctx->tail < Z_EROFS_COMPR_QUEUE_SZ) {
			ictx->pclustersize = roundup(compressedsize, blksz);
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
		if (ctx->membuf) {
			erofs_dbg("Writing %u uncompressed data of %u blocks of %s",
				  e->length, e->compressedblks, inode->i_srcpath);

			memcpy(ctx->membuf + ctx->memoff, dst - padding,
			       e->compressedblks * blksz);
			ctx->memoff += e->compressedblks * blksz;
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
	if (ctx->blkaddr != EROFS_NULL_ADDR)
		ctx->blkaddr += e->compressedblks;
	if (!may_inline && !may_packing && !is_packed_inode)
		(void)z_erofs_dedupe_insert(e, ctx->queue + ctx->head);
	ctx->head += e->length;
	return 0;

fix_dedupedfrag:
	DBG_BUGON(!inode->fragment_size);
	ctx->remaining += inode->fragment_size;
	ictx->fix_dedupedfrag = true;
	return 1;
}

static int z_erofs_compress_one(struct z_erofs_compress_sctx *ctx)
{
	struct z_erofs_compress_ictx *ictx = ctx->ictx;
	unsigned int len = ctx->tail - ctx->head;
	struct z_erofs_extent_item *ei;

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
		if (ictx->fix_dedupedfrag && !ictx->fragemitted &&
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

int z_erofs_compress_segment(struct z_erofs_compress_sctx *ctx,
			     u64 offset, erofs_blk_t blkaddr)
{
	struct z_erofs_compress_ictx *ictx = ctx->ictx;
	int fd = ictx->fd;

	ctx->blkaddr = blkaddr;
	while (ctx->remaining) {
		const u64 rx = min_t(u64, ctx->remaining,
				     Z_EROFS_COMPR_QUEUE_SZ - ctx->tail);
		int ret;

		ret = (offset == -1 ?
			read(fd, ctx->queue + ctx->tail, rx) :
			pread(fd, ctx->queue + ctx->tail, rx,
			      ictx->fpos + offset));
		if (ret != rx)
			return -errno;

		ctx->remaining -= rx;
		ctx->tail += rx;
		if (offset != -1)
			offset += rx;

		ret = z_erofs_compress_one(ctx);
		if (ret)
			return ret;
	}
	DBG_BUGON(ctx->head != ctx->tail);

	if (ctx->pivot) {
		z_erofs_commit_extent(ctx, ctx->pivot);
		ctx->pivot = NULL;
	}
	return 0;
}

int z_erofs_finalize_compression(struct z_erofs_compress_ictx *ictx,
				 struct erofs_buffer_head *bh,
				 erofs_blk_t blkaddr,
				 erofs_blk_t compressed_blocks)
{
	struct erofs_inode *inode = ictx->inode;
	struct erofs_sb_info *sbi = inode->sbi;
	u8 *compressmeta = ictx->metacur - Z_EROFS_LEGACY_MAP_HEADER_SIZE;
	unsigned int legacymetasize;
	int ret = 0;

	/* fall back to no compression mode */
	DBG_BUGON(compressed_blocks < !!inode->idata_size);
	compressed_blocks -= !!inode->idata_size;

	z_erofs_write_indexes(ictx);
	legacymetasize = ictx->metacur - compressmeta;
	/* estimate if data compression saves space or not */
	if (!inode->fragment_size &&
	    compressed_blocks * erofs_blksiz(sbi) + inode->idata_size +
	    legacymetasize >= inode->i_size) {
		z_erofs_dedupe_commit(true);

		if (inode->idata) {
			free(inode->idata);
			inode->idata = NULL;
		}
		erofs_bdrop(bh, true); /* revoke buffer */
		free(compressmeta);
		inode->compressmeta = NULL;

		return -ENOSPC;
	}
	z_erofs_dedupe_commit(false);
	z_erofs_write_mapheader(inode, compressmeta);

	if (!ictx->fragemitted)
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
}

#ifdef EROFS_MT_ENABLED
void *z_erofs_mt_wq_tls_alloc(struct erofs_workqueue *wq, void *ptr)
{
	struct erofs_compress_wq_tls *tls;

	tls = calloc(1, sizeof(*tls));
	if (!tls)
		return NULL;

	tls->queue = malloc(Z_EROFS_COMPR_QUEUE_SZ);
	if (!tls->queue)
		goto err_free_priv;

	tls->destbuf = calloc(1, EROFS_CONFIG_COMPR_MAX_SZ +
			      EROFS_MAX_BLOCK_SIZE);
	if (!tls->destbuf)
		goto err_free_queue;

	tls->ccfg = calloc(EROFS_MAX_COMPR_CFGS, sizeof(*tls->ccfg));
	if (!tls->ccfg)
		goto err_free_destbuf;
	return tls;

err_free_destbuf:
	free(tls->destbuf);
err_free_queue:
	free(tls->queue);
err_free_priv:
	free(tls);
	return NULL;
}

int z_erofs_mt_wq_tls_init_compr(struct erofs_sb_info *sbi,
				 struct erofs_compress_wq_tls *tls,
				 unsigned int alg_id, char *alg_name,
				 unsigned int comp_level,
				 unsigned int dict_size)
{
	struct erofs_compress_cfg *lc = &tls->ccfg[alg_id];
	int ret;

	if (likely(lc->enable))
		return 0;

	ret = erofs_compressor_init(sbi, &lc->handle, alg_name,
				    comp_level, dict_size);
	if (ret)
		return ret;
	lc->algorithmtype = alg_id;
	lc->enable = true;
	return 0;
}

void *z_erofs_mt_wq_tls_free(struct erofs_workqueue *wq, void *priv)
{
	struct erofs_compress_wq_tls *tls = priv;
	int i;

	for (i = 0; i < EROFS_MAX_COMPR_CFGS; i++)
		if (tls->ccfg[i].enable)
			erofs_compressor_exit(&tls->ccfg[i].handle);

	free(tls->ccfg);
	free(tls->destbuf);
	free(tls->queue);
	free(tls);
	return NULL;
}

void z_erofs_mt_workfn(struct erofs_work *work, void *tlsp)
{
	struct erofs_compress_work *cwork = (struct erofs_compress_work *)work;
	struct erofs_compress_wq_tls *tls = tlsp;
	struct z_erofs_compress_sctx *sctx = &cwork->ctx;
	struct z_erofs_mt_file *mtfile_desc = cwork->mtfile_desc;
	struct erofs_sb_info *sbi = sctx->ictx->inode->sbi;
	int ret = 0;

	ret = z_erofs_mt_wq_tls_init_compr(sbi, tls, cwork->alg_id,
					   cwork->alg_name, cwork->comp_level,
					   cwork->dict_size);
	if (ret)
		goto out;

	sctx->queue = tls->queue;
	sctx->destbuf = tls->destbuf;
	sctx->chandle = &tls->ccfg[cwork->alg_id].handle;

	sctx->membuf = malloc(round_up(sctx->remaining, erofs_blksiz(sbi)));
	if (!sctx->membuf) {
		ret = -ENOMEM;
		goto out;
	}
	sctx->memoff = 0;

	ret = z_erofs_compress_segment(sctx, sctx->seg_idx * cfg.c_segment_size,
				       EROFS_NULL_ADDR);

out:
	cwork->errcode = ret;
	pthread_mutex_lock(&mtfile_desc->mutex);
	++mtfile_desc->nfini;
	pthread_cond_signal(&mtfile_desc->cond);
	pthread_mutex_unlock(&mtfile_desc->mutex);
}

int z_erofs_merge_segment(struct z_erofs_compress_ictx *ictx,
			  struct z_erofs_compress_sctx *sctx)
{
	struct z_erofs_extent_item *ei, *n;
	struct erofs_sb_info *sbi = ictx->inode->sbi;
	erofs_blk_t blkoff = 0;
	int ret = 0, ret2;

	list_for_each_entry_safe(ei, n, &sctx->extents, list) {
		list_del(&ei->list);
		list_add_tail(&ei->list, &ictx->extents);

		if (ei->e.blkaddr != EROFS_NULL_ADDR)	/* deduped extents */
			continue;

		ei->e.blkaddr = sctx->blkaddr;
		sctx->blkaddr += ei->e.compressedblks;

		ret2 = blk_write(sbi, sctx->membuf + blkoff * erofs_blksiz(sbi),
				 ei->e.blkaddr, ei->e.compressedblks);
		blkoff += ei->e.compressedblks;
		if (ret2) {
			ret = ret2;
			continue;
		}
	}
	free(sctx->membuf);
	return ret;
}

int z_erofs_mt_compress(struct z_erofs_compress_ictx *ictx,
			struct erofs_compress_cfg *ccfg)
{
	struct erofs_compress_work *cur, *head = NULL, **last = &head;
	struct erofs_inode *inode = ictx->inode;
	struct z_erofs_mt_file *mtfile_desc = NULL;
	int nsegs = DIV_ROUND_UP(inode->i_size, cfg.c_segment_size);
	int i;

	pthread_mutex_lock(&z_erofs_mt_ctrl.file_mutex);
	if (z_erofs_mt_ctrl.file_idle) {
		mtfile_desc = z_erofs_mt_ctrl.file_idle;
		z_erofs_mt_ctrl.file_idle = mtfile_desc->next;
		mtfile_desc->next = NULL;
	}
	pthread_mutex_unlock(&z_erofs_mt_ctrl.file_mutex);
	if (!mtfile_desc) {
		mtfile_desc = calloc(1, sizeof(*mtfile_desc));
		if (!mtfile_desc) {
			return -ENOMEM;
		}
	}
	inode->mt_desc = mtfile_desc;

	mtfile_desc->fd = ictx->fd;
	mtfile_desc->total = nsegs;
	mtfile_desc->nfini = 0;
	pthread_mutex_init(&mtfile_desc->mutex, NULL);
	pthread_cond_init(&mtfile_desc->cond, NULL);

	for (i = 0; i < nsegs; i++) {
		cur = NULL;

		pthread_mutex_lock(&z_erofs_mt_ctrl.work_mutex);
		if (z_erofs_mt_ctrl.work_idle) {
			cur = z_erofs_mt_ctrl.work_idle;
			z_erofs_mt_ctrl.work_idle = cur->next;
			cur->next = NULL;
		}
		pthread_mutex_unlock(&z_erofs_mt_ctrl.work_mutex);
		if (!cur) {
			cur = calloc(1, sizeof(*cur));
			if (!cur) {
				while (head) {
					cur = head;
					head = cur->next;
					free(cur);
				}
				free(mtfile_desc);
				return -ENOMEM;
			}
		}

		if (i == 0)
			mtfile_desc->head = cur;
		*last = cur;
		last = &cur->next;

		cur->ctx = (struct z_erofs_compress_sctx) {
			.ictx = ictx,
			.seg_num = nsegs,
			.seg_idx = i,
			.pivot = &dummy_pivot,
		};
		init_list_head(&cur->ctx.extents);

		if (i == nsegs - 1)
			cur->ctx.remaining = inode->i_size -
					      inode->fragment_size -
					      i * cfg.c_segment_size;
		else
			cur->ctx.remaining = cfg.c_segment_size;

		cur->alg_id = ccfg->handle.alg->id;
		cur->alg_name = ccfg->handle.alg->name;
		cur->comp_level = ccfg->handle.compression_level;
		cur->dict_size = ccfg->handle.dict_size;

		cur->mtfile_desc = mtfile_desc;
		cur->work.fn = z_erofs_mt_workfn;
		erofs_queue_work(&z_erofs_mt_ctrl.wq, &cur->work);
	}

	return 0;
}

int z_erofs_mt_reap(struct z_erofs_mt_file *desc) {
	struct erofs_buffer_head *bh = NULL;
	struct erofs_compress_work *cur = desc->head, *tmp;
	struct z_erofs_compress_ictx *ictx = cur->ctx.ictx;
	erofs_blk_t blkaddr, compressed_blocks = 0;
	int ret = 0;

	bh = erofs_balloc(DATA, 0, 0, 0);
	if (IS_ERR(bh)) {
		ret = PTR_ERR(bh);
		goto out;
	}
	blkaddr = erofs_mapbh(bh->block);

	while (cur) {
		if (cur->errcode) {
			ret = cur->errcode;
		} else {
			int ret2;

			cur->ctx.blkaddr = blkaddr;
			ret2 = z_erofs_merge_segment(ictx, &cur->ctx);
			if (ret2)
				ret = ret2;

			compressed_blocks += cur->ctx.blkaddr - blkaddr;
			blkaddr = cur->ctx.blkaddr;
		}

		tmp = cur->next;
		pthread_mutex_lock(&z_erofs_mt_ctrl.work_mutex);
		cur->next = z_erofs_mt_ctrl.work_idle;
		z_erofs_mt_ctrl.work_idle = cur;
		pthread_mutex_unlock(&z_erofs_mt_ctrl.work_mutex);
		cur = tmp;
	}
	if (ret)
		goto out;

	ret = z_erofs_finalize_compression(
		ictx, bh, blkaddr - compressed_blocks, compressed_blocks);

out:
	free(ictx);
	pthread_mutex_lock(&z_erofs_mt_ctrl.file_mutex);
	desc->next = z_erofs_mt_ctrl.file_idle;
	z_erofs_mt_ctrl.file_idle = desc;
	pthread_mutex_unlock(&z_erofs_mt_ctrl.file_mutex);

	return ret;
}
#endif

int erofs_write_compressed_file(struct erofs_inode *inode, int fd, u64 fpos)
{
	static u8 g_queue[Z_EROFS_COMPR_QUEUE_SZ];
	struct erofs_buffer_head *bh;
	static struct z_erofs_compress_ictx ctx;
	static struct z_erofs_compress_sctx sctx;
	struct erofs_compress_cfg *ccfg;
	erofs_blk_t blkaddr, compressed_blocks = 0;
	int ret;
	struct erofs_sb_info *sbi = inode->sbi;
	u8 *compressmeta = malloc(BLK_ROUND_UP(sbi, inode->i_size) *
				  sizeof(struct z_erofs_lcluster_index) +
				  Z_EROFS_LEGACY_MAP_HEADER_SIZE);

	if (!compressmeta)
		return -ENOMEM;

	if (!z_erofs_mt_enabled) {
		/* allocate main data buffer */
		bh = erofs_balloc(DATA, 0, 0, 0);
		if (IS_ERR(bh)) {
			ret = PTR_ERR(bh);
			goto err_free_meta;
		}
		blkaddr = erofs_mapbh(bh->block); /* start_blkaddr */
	} else {
		bh = NULL;
		blkaddr = EROFS_NULL_ADDR;
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
		ret = z_erofs_fragments_dedupe(inode, fd, &ctx.tof_chksum);
		if (ret < 0)
			goto err_bdrop;
	}

	ctx.inode = inode;
	ctx.pclustersize = z_erofs_get_max_pclustersize(inode);
	ctx.fd = fd;
	ctx.fpos = fpos;
	ctx.metacur = compressmeta + Z_EROFS_LEGACY_MAP_HEADER_SIZE;
	init_list_head(&ctx.extents);
	ctx.fix_dedupedfrag = false;
	ctx.fragemitted = false;
	sctx = (struct z_erofs_compress_sctx) { .ictx = &ctx, };
	init_list_head(&sctx.extents);

	if (cfg.c_all_fragments && !erofs_is_packed_inode(inode) &&
	    !inode->fragment_size) {
		ret = z_erofs_pack_file_from_fd(inode, fd, ctx.tof_chksum);
		if (ret)
			goto err_free_idata;
#ifdef EROFS_MT_ENABLED
	} else if (z_erofs_mt_enabled) {
		struct z_erofs_compress_ictx *l_ictx;

		l_ictx = malloc(sizeof(*l_ictx));
		if (!l_ictx) {
			ret = -ENOMEM;
			goto err_free_idata;
		}

		memcpy(l_ictx, &ctx, sizeof(*l_ictx));
		init_list_head(&l_ictx->extents);

		ret = z_erofs_mt_compress(l_ictx, ccfg);
		if (ret) {
			free(l_ictx);
			goto err_free_idata;
		}
		return 0;
#endif
	} else {
		sctx.queue = g_queue;
		sctx.destbuf = NULL;
		sctx.chandle = &ccfg->handle;
		sctx.remaining = inode->i_size - inode->fragment_size;
		sctx.seg_num = 1;
		sctx.seg_idx = 0;
		sctx.pivot = &dummy_pivot;

		ret = z_erofs_compress_segment(&sctx, -1, blkaddr);
		if (ret)
			goto err_free_idata;
		compressed_blocks = sctx.blkaddr - blkaddr;
	}

	/* generate an extent for the deduplicated fragment */
	if (inode->fragment_size && !ctx.fragemitted) {
		struct z_erofs_extent_item *ei;

		ei = malloc(sizeof(*ei));
		if (!ei) {
			ret = -ENOMEM;
			goto err_free_idata;
		}

		ei->e = (struct z_erofs_inmem_extent) {
			.length = inode->fragment_size,
			.compressedblks = 0,
			.raw = false,
			.partial = false,
			.blkaddr = sctx.blkaddr,
		};
		init_list_head(&ei->list);
		z_erofs_commit_extent(&sctx, ei);
	}
	z_erofs_fragments_commit(inode);
	list_splice_tail(&sctx.extents, &ctx.extents);

	return z_erofs_finalize_compression(&ctx, bh, blkaddr,
					    compressed_blocks);

err_free_idata:
	if (inode->idata) {
		free(inode->idata);
		inode->idata = NULL;
	}
err_bdrop:
	if (bh)
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

	z_erofs_mt_enabled = false;
#ifdef EROFS_MT_ENABLED
	if (cfg.c_mt_workers > 1) {
		pthread_mutex_init(&z_erofs_mt_ctrl.file_mutex, NULL);
		pthread_mutex_init(&z_erofs_mt_ctrl.work_mutex, NULL);
		ret = erofs_alloc_workqueue(&z_erofs_mt_ctrl.wq,
					    cfg.c_mt_workers,
					    cfg.c_mt_workers << 2,
					    z_erofs_mt_wq_tls_alloc,
					    z_erofs_mt_wq_tls_free);
		z_erofs_mt_enabled = !ret;
	}
#endif
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
		while (z_erofs_mt_ctrl.work_idle) {
			struct erofs_compress_work *tmp =
				z_erofs_mt_ctrl.work_idle->next;
			free(z_erofs_mt_ctrl.work_idle);
			z_erofs_mt_ctrl.work_idle = tmp;
		}
		while (z_erofs_mt_ctrl.file_idle) {
			struct z_erofs_mt_file *tmp =
				z_erofs_mt_ctrl.file_idle->next;
			free(z_erofs_mt_ctrl.file_idle);
			z_erofs_mt_ctrl.file_idle = tmp;
		}
#endif
	}
	return 0;
}
