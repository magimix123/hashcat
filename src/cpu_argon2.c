/**
 * Author......: See docs/credits.txt
 * License.....: MIT
 */

#define IS_GENERIC

#include "common.h"
#include "types.h"
#include "bitops.h"
#include "inc_hash_constants.h"
#include "inc_hash_functions.cl"
#include "cpu_argon2.h"


static argon2_ctx_t *argon2_ctx = NULL;

void argon2_initialize(argon2_t *argon2_params)
{
  if (argon2_ctx == NULL) argon2_ctx = (argon2_ctx_t *) hcmalloc(sizeof(argon2_ctx_t));


}

void argon2_destory()
{
  if (argon2_ctx != NULL) free (argon2_ctx);
}
