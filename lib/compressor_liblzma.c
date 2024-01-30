// SPDX-License-Identifier: GPL-2.0+ OR Apache-2.0
/*
 * Copyright (C) 2021 Gao Xiang <xiang@kernel.org>
 */
#include <stdlib.h>
#include "config.h"
#ifdef HAVE_LIBLZMA
#include <lzma.h>
#include "erofs/config.h"
#include "erofs/print.h"
#include "erofs/internal.h"
#include "compressor.h"

struct erofs_liblzma_context {
	lzma_options_lzma opt;
	lzma_stream strm;
};

static int erofs_liblzma_compress_destsize(const struct erofs_compress *c,
					   const void *src, unsigned int *srcsize,
					   void *dst, unsigned int dstsize)
{
	struct erofs_liblzma_context *ctx = c->private_data;
	lzma_stream *strm = &ctx->strm;

	lzma_ret ret = lzma_microlzma_encoder(strm, &ctx->opt);
	if (ret != LZMA_OK)
		return -EFAULT;

	strm->next_in = src;
	strm->avail_in = *srcsize;
	strm->next_out = dst;
	strm->avail_out = dstsize;

	ret = lzma_code(strm, LZMA_FINISH);
	if (ret != LZMA_STREAM_END)
		return -EBADMSG;

	*srcsize = strm->total_in;
	return strm->total_out;
}

static int erofs_compressor_liblzma_exit(struct erofs_compress *c)
{
	struct erofs_liblzma_context *ctx = c->private_data;

	if (!ctx)
		return -EINVAL;

	lzma_end(&ctx->strm);
	free(ctx);
	return 0;
}

static int erofs_compressor_liblzma_setlevel(struct erofs_compress *c,
					     int compression_level)
{
	if (compression_level < 0)
		compression_level = erofs_compressor_lzma.default_level;

	if (compression_level > erofs_compressor_lzma.best_level) {
		erofs_err("invalid compression level %d", compression_level);
		return -EINVAL;
	}
	c->compression_level = compression_level;
	return 0;
}

static int erofs_compressor_liblzma_setdictsize(struct erofs_compress *c,
						u32 dict_size)
{
	if (!dict_size)
		dict_size = erofs_compressor_lzma.default_dictsize;

	if (dict_size > erofs_compressor_lzma.max_dictsize ||
	    dict_size < 4096) {
		erofs_err("invalid dictionary size %u", dict_size);
		return -EINVAL;
	}
	c->dict_size = dict_size;
	return 0;
}

static int erofs_compressor_liblzma_init(struct erofs_compress *c,
					 bool print_warning)
{
	struct erofs_liblzma_context *ctx;
	u32 preset;

	ctx = malloc(sizeof(*ctx));
	if (!ctx)
		return -ENOMEM;
	ctx->strm = (lzma_stream)LZMA_STREAM_INIT;

	if (c->compression_level < 0)
		preset = LZMA_PRESET_DEFAULT;
	else if (c->compression_level >= 100)
		preset = (c->compression_level - 100) | LZMA_PRESET_EXTREME;
	else
		preset = c->compression_level;

	if (lzma_lzma_preset(&ctx->opt, preset))
		return -EINVAL;
	ctx->opt.dict_size = c->dict_size;

	c->private_data = ctx;
	if (print_warning) {
		erofs_warn("EXPERIMENTAL MicroLZMA feature in use. Use at your own risk!");
		erofs_warn("Note that it may take more time since the compressor is still single-threaded for now.");
	}
	return 0;
}

const struct erofs_compressor erofs_compressor_lzma = {
	.default_level = LZMA_PRESET_DEFAULT,
	.best_level = 109,
	.default_dictsize = Z_EROFS_LZMA_MAX_DICT_SIZE,
	.max_dictsize = Z_EROFS_LZMA_MAX_DICT_SIZE,
	.init = erofs_compressor_liblzma_init,
	.exit = erofs_compressor_liblzma_exit,
	.setlevel = erofs_compressor_liblzma_setlevel,
	.setdictsize = erofs_compressor_liblzma_setdictsize,
	.compress_destsize = erofs_liblzma_compress_destsize,
};
#endif
