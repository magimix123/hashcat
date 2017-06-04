/**
 * Author......: See docs/credits.txt
 * License.....: MIT
 */

#define NEW_SIMD_CODE

#include "inc_vendor.cl"
#include "inc_hash_constants.h"
#include "inc_hash_functions.cl"
#include "inc_types.cl"
#include "inc_common.cl"
#include "inc_rp.h"
#include "inc_rp.cl"
#include "inc_simd.cl"

#define GRAIN128A_IVSIZE    96
#define GRAIN128A_KEYSIZE  128

void grain_keystream(u8x LFSR[128], u8x NFSR[128], u8x outbit)
{
  u8x NBit, LBit;

  /* Calculate feedback and output bits */

  outbit = NFSR[2]                                     \
         ^ NFSR[15]                                    \
         ^ NFSR[36]                                    \
         ^ NFSR[45]                                    \
         ^ NFSR[64]                                    \
         ^ NFSR[73]                                    \
         ^ NFSR[89]                                    \
         ^ LFSR[93]                                    \
         ^ (NFSR[12] & LFSR[ 8])                       \
         ^ (LFSR[13] & LFSR[20])                       \
         ^ (NFSR[95] & LFSR[42])                       \
         ^ (LFSR[60] & LFSR[79])                       \
         ^ (NFSR[12] & NFSR[95] & LFSR[94]);    

  NBit   = LFSR[0]                                            \
         ^ NFSR[0]                                            \
         ^ NFSR[26]                                           \
         ^ NFSR[56]                                           \
         ^ NFSR[91]                                           \
         ^ NFSR[96]                                           \
         ^ (NFSR[3]  & NFSR[67])                              \
         ^ (NFSR[11] & NFSR[13])                              \
         ^ (NFSR[17] & NFSR[18])                              \
         ^ (NFSR[27] & NFSR[59])                              \
         ^ (NFSR[40] & NFSR[48])                              \
         ^ (NFSR[61] & NFSR[65])                              \
         ^ (NFSR[68] & NFSR[84])                              \
         ^ (NFSR[88] & NFSR[92] & NFSR[93] & NFSR[95])        \
         ^ (NFSR[22] & NFSR[24] & NFSR[25])                   \
         ^ (NFSR[70] & NFSR[78] & NFSR[82]);

  LBit   = LFSR[ 0]                         \
         ^ LFSR[ 7]                         \
         ^ LFSR[38]                         \
         ^ LFSR[70]                         \
         ^ LFSR[81]                         \
         ^ LFSR[96];

  /* Update registers */

  for (u8 i = 1; i < GRAIN128A_KEYSIZE; i++) 
  {
    NFSR[i - 1] = NFSR[i];
    LFSR[i - 1] = LFSR[i];
  }

  NFSR[GRAIN128A_KEYSIZE - 1] = NBit;
  LFSR[GRAIN128A_KEYSIZE - 1] = LBit;
}

void grain128a_transform (const u32x key[4], const u32 iv[3], const u32 plain[2], const u32x digest[8])
{
  /**
  * Initialize registers
  */

  u8x outbit;
  u8x LFSR[128];
  u8x NFSR[128];

  for (u8 i = 0; i< GRAIN128A_IVSIZE / 8; i++) 
  {
    for (u8 j = 0; j < 8; j++) 
    {
      NFSR[i * 8 + j] = ((key[i] >> j) & 1 );
      LFSR[i * 8 + j] = ((iv[i]  >> j) & 1 );
    }
  }

  for (u8 i = GRAIN128A_IVSIZE / 8; i < GRAIN128A_KEYSIZE / 8; i++) 
  {
    for (u8 j = 0; j < 8; j++) 
    {
      NFSR[i * 8 + j] = ((key[i] >> j ) & 1);
      LFSR[i * 8 + j] = 1;
    }
  }

  LFSR[127] = 0; // padding req

  /* do initial clockings */

  for (u8 i = 0; i < 256; i++) 
  {
    grain_keystream(LFSR, NFSR, outbit);

    LFSR[127] ^= outbit;
    NFSR[127] ^= outbit;             
  }

  u8 k = 0;

  for (u8 i = 0; i < 8; i++) 
  {
    k = 0;
    for (u8 j = 0; j < 8; j++) 
    {
      grain_keystream(LFSR, NFSR, outbit);
      k |= outbit << j;
    }

    //digest[i] = plain[i] ^ k;
  }
}

__kernel void m15600_m04 (__global pw_t *pws, __global const kernel_rule_t *rules_buf, __global const comb_t *combs_buf, __global const bf_t *bfs_buf, __global void *tmps, __global void *hooks, __global const u32 *bitmaps_buf_s1_a, __global const u32 *bitmaps_buf_s1_b, __global const u32 *bitmaps_buf_s1_c, __global const u32 *bitmaps_buf_s1_d, __global const u32 *bitmaps_buf_s2_a, __global const u32 *bitmaps_buf_s2_b, __global const u32 *bitmaps_buf_s2_c, __global const u32 *bitmaps_buf_s2_d, __global plain_t *plains_buf, __global const digest_t *digests_buf, __global u32 *hashes_shown, __global const salt_t *salt_bufs, __global const grain128a_t *esalt_bufs, __global u32 *d_return_buf, __global u32 *d_scryptV0_buf, __global u32 *d_scryptV1_buf, __global u32 *d_scryptV2_buf, __global u32 *d_scryptV3_buf, const u32 bitmap_mask, const u32 bitmap_shift1, const u32 bitmap_shift2, const u32 salt_pos, const u32 loop_pos, const u32 loop_cnt, const u32 il_cnt, const u32 digests_cnt, const u32 digests_offset, const u32 combs_mode, const u32 gid_max)
{ 
  /**
   * modifier
   */

  const u32 gid = get_global_id (0);
  const u32 lid = get_local_id (0);

  u32 pw_buf0[4];
  u32 pw_buf1[4];

  pw_buf0[0] = pws[gid].i[0];
  pw_buf0[1] = pws[gid].i[1];
  pw_buf0[2] = pws[gid].i[2];
  pw_buf0[3] = pws[gid].i[3];
  pw_buf1[0] = pws[gid].i[4];
  pw_buf1[1] = pws[gid].i[5];
  pw_buf1[2] = pws[gid].i[6];
  pw_buf1[3] = pws[gid].i[7];

  const u32 pw_len = pws[gid].pw_len;

  /**
   * Assign global to private       
   */

  u32 iv[3]       = { 0 };
  u32 plain[2]    = { 0 };

  iv[0] = esalt_bufs[digests_offset].iv[0];
  iv[1] = esalt_bufs[digests_offset].iv[1];
  iv[2] = esalt_bufs[digests_offset].iv[2];

  plain[0] = esalt_bufs[digests_offset].plain[0];
  plain[1] = esalt_bufs[digests_offset].plain[1];

  /**
   * loop
   */

  for (u32 il_pos = 0; il_pos < il_cnt; il_pos += VECT_SIZE)
  {
    u32x w0[4] = { 0 };
    u32x w1[4] = { 0 };

    const u32x out_len = apply_rules_vect(pw_buf0, pw_buf1, pw_len, rules_buf, il_pos, w0, w1);

    u8x  key[16]   = { 0 };
    u32x digest[4] = { 0 };

    key[ 0] = (const u8) w0[0];

    grain128a_transform (w0, iv, plain, digest);

    const u32x r0 = digest[0];
    const u32x r1 = digest[1];
    const u32x r2 = digest[2];
    const u32x r3 = digest[3];

    COMPARE_M_SIMD(r0, r1, r2, r3);
  } 
}

__kernel void m15600_m08 (__global pw_t *pws, __global const kernel_rule_t *rules_buf, __global const comb_t *combs_buf, __global const bf_t *bfs_buf, __global void *tmps, __global void *hooks, __global const u32 *bitmaps_buf_s1_a, __global const u32 *bitmaps_buf_s1_b, __global const u32 *bitmaps_buf_s1_c, __global const u32 *bitmaps_buf_s1_d, __global const u32 *bitmaps_buf_s2_a, __global const u32 *bitmaps_buf_s2_b, __global const u32 *bitmaps_buf_s2_c, __global const u32 *bitmaps_buf_s2_d, __global plain_t *plains_buf, __global const digest_t *digests_buf, __global u32 *hashes_shown, __global const salt_t *salt_bufs, __global const grain128a_t *esalt_bufs, __global u32 *d_return_buf, __global u32 *d_scryptV0_buf, __global u32 *d_scryptV1_buf, __global u32 *d_scryptV2_buf, __global u32 *d_scryptV3_buf, const u32 bitmap_mask, const u32 bitmap_shift1, const u32 bitmap_shift2, const u32 salt_pos, const u32 loop_pos, const u32 loop_cnt, const u32 il_cnt, const u32 digests_cnt, const u32 digests_offset, const u32 combs_mode, const u32 gid_max)
{
}

__kernel void m15600_m16 (__global pw_t *pws, __global const kernel_rule_t *rules_buf, __global const comb_t *combs_buf, __global const bf_t *bfs_buf, __global void *tmps, __global void *hooks, __global const u32 *bitmaps_buf_s1_a, __global const u32 *bitmaps_buf_s1_b, __global const u32 *bitmaps_buf_s1_c, __global const u32 *bitmaps_buf_s1_d, __global const u32 *bitmaps_buf_s2_a, __global const u32 *bitmaps_buf_s2_b, __global const u32 *bitmaps_buf_s2_c, __global const u32 *bitmaps_buf_s2_d, __global plain_t *plains_buf, __global const digest_t *digests_buf, __global u32 *hashes_shown, __global const salt_t *salt_bufs, __global const grain128a_t *esalt_bufs, __global u32 *d_return_buf, __global u32 *d_scryptV0_buf, __global u32 *d_scryptV1_buf, __global u32 *d_scryptV2_buf, __global u32 *d_scryptV3_buf, const u32 bitmap_mask, const u32 bitmap_shift1, const u32 bitmap_shift2, const u32 salt_pos, const u32 loop_pos, const u32 loop_cnt, const u32 il_cnt, const u32 digests_cnt, const u32 digests_offset, const u32 combs_mode, const u32 gid_max)
{
}

__kernel void m15600_s04 (__global pw_t *pws, __global const kernel_rule_t *rules_buf, __global const comb_t *combs_buf, __global const bf_t *bfs_buf, __global void *tmps, __global void *hooks, __global const u32 *bitmaps_buf_s1_a, __global const u32 *bitmaps_buf_s1_b, __global const u32 *bitmaps_buf_s1_c, __global const u32 *bitmaps_buf_s1_d, __global const u32 *bitmaps_buf_s2_a, __global const u32 *bitmaps_buf_s2_b, __global const u32 *bitmaps_buf_s2_c, __global const u32 *bitmaps_buf_s2_d, __global plain_t *plains_buf, __global const digest_t *digests_buf, __global u32 *hashes_shown, __global const salt_t *salt_bufs, __global const grain128a_t *esalt_bufs, __global u32 *d_return_buf, __global u32 *d_scryptV0_buf, __global u32 *d_scryptV1_buf, __global u32 *d_scryptV2_buf, __global u32 *d_scryptV3_buf, const u32 bitmap_mask, const u32 bitmap_shift1, const u32 bitmap_shift2, const u32 salt_pos, const u32 loop_pos, const u32 loop_cnt, const u32 il_cnt, const u32 digests_cnt, const u32 digests_offset, const u32 combs_mode, const u32 gid_max)
{ 
  /**
   * modifier
   */

  const u32 lid = get_local_id (0);

  const u32 gid = get_global_id (0);

  if (gid >= gid_max) return;

  u32 pw_buf0[4];
  u32 pw_buf1[4];

  pw_buf0[0] = pws[gid].i[0];
  pw_buf0[1] = pws[gid].i[1];
  pw_buf0[2] = pws[gid].i[2];
  pw_buf0[3] = pws[gid].i[3];
  pw_buf1[0] = pws[gid].i[4];
  pw_buf1[1] = pws[gid].i[5];
  pw_buf1[2] = pws[gid].i[6];
  pw_buf1[3] = pws[gid].i[7];

  const u32 pw_len = pws[gid].pw_len;

  /**
   * Assign global to private
   */

  u32 iv[3]       = { 0 };
  u32 plain[2]    = { 0 };

  iv[0] = esalt_bufs[digests_offset].iv[0];
  iv[1] = esalt_bufs[digests_offset].iv[1];
  iv[2] = esalt_bufs[digests_offset].iv[2];

  plain[0] = esalt_bufs[digests_offset].plain[0];
  plain[1] = esalt_bufs[digests_offset].plain[1];
  
  /**
   * digest
   */

  const u32 search[4] =
  {
    digests_buf[digests_offset].digest_buf[DGST_R0],
    digests_buf[digests_offset].digest_buf[DGST_R1],
    digests_buf[digests_offset].digest_buf[DGST_R2],
    digests_buf[digests_offset].digest_buf[DGST_R3]
  };
 
  /**
   * loop
   */

  for (u32 il_pos = 0; il_pos < il_cnt; il_pos += VECT_SIZE)
  {
    u32x w0[4] = { 0 };
    u32x w1[4] = { 0 };
    
    const u32x out_len = apply_rules_vect(pw_buf0, pw_buf1, pw_len, rules_buf, il_pos, w0, w1);

    u32x digest[4] = { 0 };

//    grain128a_transform (w0, iv, plain, digest);

    const u32x r0 = digest[0];
    const u32x r1 = digest[1];
    const u32x r2 = digest[2];
    const u32x r3 = digest[3];

    COMPARE_S_SIMD(r0, r1, r2, r3);
  }  
}

__kernel void m15600_s08 (__global pw_t *pws, __global const kernel_rule_t *rules_buf, __global const comb_t *combs_buf, __global const bf_t *bfs_buf, __global void *tmps, __global void *hooks, __global const u32 *bitmaps_buf_s1_a, __global const u32 *bitmaps_buf_s1_b, __global const u32 *bitmaps_buf_s1_c, __global const u32 *bitmaps_buf_s1_d, __global const u32 *bitmaps_buf_s2_a, __global const u32 *bitmaps_buf_s2_b, __global const u32 *bitmaps_buf_s2_c, __global const u32 *bitmaps_buf_s2_d, __global plain_t *plains_buf, __global const digest_t *digests_buf, __global u32 *hashes_shown, __global const salt_t *salt_bufs, __global const grain128a_t *esalt_bufs, __global u32 *d_return_buf, __global u32 *d_scryptV0_buf, __global u32 *d_scryptV1_buf, __global u32 *d_scryptV2_buf, __global u32 *d_scryptV3_buf, const u32 bitmap_mask, const u32 bitmap_shift1, const u32 bitmap_shift2, const u32 salt_pos, const u32 loop_pos, const u32 loop_cnt, const u32 il_cnt, const u32 digests_cnt, const u32 digests_offset, const u32 combs_mode, const u32 gid_max)
{
}

__kernel void m15600_s16 (__global pw_t *pws, __global const kernel_rule_t *rules_buf, __global const comb_t *combs_buf, __global const bf_t *bfs_buf, __global void *tmps, __global void *hooks, __global const u32 *bitmaps_buf_s1_a, __global const u32 *bitmaps_buf_s1_b, __global const u32 *bitmaps_buf_s1_c, __global const u32 *bitmaps_buf_s1_d, __global const u32 *bitmaps_buf_s2_a, __global const u32 *bitmaps_buf_s2_b, __global const u32 *bitmaps_buf_s2_c, __global const u32 *bitmaps_buf_s2_d, __global plain_t *plains_buf, __global const digest_t *digests_buf, __global u32 *hashes_shown, __global const salt_t *salt_bufs, __global const grain128a_t *esalt_bufs, __global u32 *d_return_buf, __global u32 *d_scryptV0_buf, __global u32 *d_scryptV1_buf, __global u32 *d_scryptV2_buf, __global u32 *d_scryptV3_buf, const u32 bitmap_mask, const u32 bitmap_shift1, const u32 bitmap_shift2, const u32 salt_pos, const u32 loop_pos, const u32 loop_cnt, const u32 il_cnt, const u32 digests_cnt, const u32 digests_offset, const u32 combs_mode, const u32 gid_max)
{
}
