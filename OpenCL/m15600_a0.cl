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

#define GRAIN128A_IVSIZE            96
#define GRAIN128A_KEYSIZE          128
#define GRAIN128A_LFSR_PAD  0xfffffffe

#define update_keystream()                                                              \
  do                                                                                    \
  {                                                                                     \
    outbit = ((NFSR[0]  >>  2) & 1)                                                     \
           ^ ((NFSR[0]  >> 15) & 1)                                                     \
           ^ ((NFSR[1]  >>  4) & 1)                                                     \
           ^ ((NFSR[1]  >> 13) & 1)                                                     \
           ^ ((NFSR[2]  >>  0) & 1)                                                     \
           ^ ((NFSR[2]  >>  9) & 1)                                                     \
           ^ ((NFSR[2]  >> 25) & 1)                                                     \
           ^ ((LFSR[2]  >> 29) & 1)                                                     \
           ^ (((NFSR[0] >> 12) & 1) & ((LFSR[0] >>  8) & 1))                            \
           ^ (((LFSR[0] >> 13) & 1) & ((LFSR[0] >> 20) & 1))                            \
           ^ (((NFSR[2] >> 31) & 1) & ((LFSR[1] >> 10) & 1))                            \
           ^ (((LFSR[1] >> 28) & 1) & ((LFSR[2] >> 15) & 1))                            \
           ^ (((NFSR[0] >> 12) & 1) & ((NFSR[2] >> 31) & 1) & ((LFSR[2] >> 30) & 1));   \
                                                                                        \
    NBit   = ((LFSR[0] >>  0) & 1)                                                                              \
           ^ ((NFSR[0] >>  0) & 1)                                                                              \
           ^ ((NFSR[0] >> 26) & 1)                                                                              \
           ^ ((NFSR[1] >> 24) & 1)                                                                              \
           ^ ((NFSR[2] >> 27) & 1)                                                                              \
           ^ ((NFSR[3] >>  0) & 1)                                                                              \
           ^ (((NFSR[0] >>  3) & 1) & ((NFSR[2] >>  3) & 1))                                                    \
           ^ (((NFSR[0] >> 11) & 1) & ((NFSR[0] >> 13) & 1))                                                    \
           ^ (((NFSR[0] >> 17) & 1) & ((NFSR[0] >> 18) & 1))                                                    \
           ^ (((NFSR[0] >> 27) & 1) & ((NFSR[1] >> 27) & 1))                                                    \
           ^ (((NFSR[1] >>  8) & 1) & ((NFSR[1] >> 16) & 1))                                                    \
           ^ (((NFSR[1] >> 29) & 1) & ((NFSR[2] >>  1) & 1))                                                    \
           ^ (((NFSR[2] >>  4) & 1) & ((NFSR[2] >> 20) & 1))                                                    \
           ^ (((NFSR[2] >> 24) & 1) & ((NFSR[2] >> 28) & 1) & ((NFSR[2] >> 29) & 1) & ((NFSR[2] >> 31) & 1))    \
           ^ (((NFSR[0] >> 22) & 1) & ((NFSR[0] >> 24) & 1) & ((NFSR[0] >> 25) & 1))                            \
           ^ (((NFSR[2] >>  6) & 1) & ((NFSR[2] >> 14) & 1) & ((NFSR[2] >> 18) & 1));                           \
                                                                                                                \
    LBit   = ((LFSR[0] >>  0) & 1)                      \
           ^ ((LFSR[0] >>  7) & 1)                      \
           ^ ((LFSR[1] >>  6) & 1)                      \
           ^ ((LFSR[2] >>  6) & 1)                      \
           ^ ((LFSR[2] >> 17) & 1)                      \
           ^ ((LFSR[3] >>  0) & 1);                     \
                                                        \
    NFSR[0] >>= 1;                                      \
    NFSR[0]  ^= ((NFSR[1] << 31) & 0x80000000);         \
    NFSR[1] >>= 1;                                      \
    NFSR[1]  ^= ((NFSR[2] << 31) & 0x80000000);         \
    NFSR[2] >>= 1;                                      \
    NFSR[2]  ^= ((NFSR[3] << 31) & 0x80000000);         \
    NFSR[3] >>= 1;                                      \
    NFSR[3]  ^= ((NBit    << 31) & 0x80000000);         \
    LFSR[0] >>= 1;                                      \
    LFSR[0]  ^= ((LFSR[1] << 31) & 0x80000000);         \
    LFSR[1] >>= 1;                                      \
    LFSR[1]  ^= ((LFSR[2] << 31) & 0x80000000);         \
    LFSR[2] >>= 1;                                      \
    LFSR[2]  ^= ((LFSR[3] << 31) & 0x80000000);         \
    LFSR[3] >>= 1;                                      \
    LFSR[3]  ^= ((LBit    << 31) & 0x80000000);         \
                                                        \
  } while (0)

void grain128a_transform (const u32x key[4], const u32 iv[3], const u32 plain[2], u8x out[8])
{
  /**
   * Initialize registers
   */
 
  u32x  outbit  = 0;

  u32x LFSR[4];
  u32x NFSR[4];
  u32x LBit;
  u32x NBit;

  LFSR[0] = iv[0];
  LFSR[1] = iv[1];
  LFSR[2] = iv[2];
  LFSR[3] = GRAIN128A_LFSR_PAD;  

  NFSR[0] = key[0];
  NFSR[1] = key[1];
  NFSR[2] = key[2];
  NFSR[3] = key[3];
  
  /* do initial clockings */

  for (u8 i = 0; i < 128; i++)
  {
    update_keystream();

    LFSR[3] ^= ((outbit << 31) & 0x80000000);
    NFSR[3] ^= ((outbit << 31) & 0x80000000);
  }

  for (u8 i = 0; i < 128; i++)
  {
    update_keystream();

    LFSR[3] ^= ((outbit << 31) & 0x80000000);
    NFSR[3] ^= ((outbit << 31) & 0x80000000);
  }

  u32 k = 0;

  for (u8 i = 0; i < 8; i++)
  {

    k = 0;
    for (u8 j = 0; j < 8; j++)
    {
      update_keystream();
      k ^= (outbit << j);
    }
    out[i] = k;
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
    u32x testkey[4] = { 0 };
    
    const u32x out_len = apply_rules_vect(pw_buf0, pw_buf1, pw_len, rules_buf, il_pos, w0, w1);

    u32x digest[4] = { 0 };
    u8x out[8]    = { 0 };

    testkey[0] = 0x00000000;
    testkey[1] = 0x00000000;
    testkey[2] = 0x00000000;
    testkey[3] = 0x00000000;

    grain128a_transform (w0, iv, plain, out);

    printf("out: %02x%02x%02x%02x%02x%02x%02x%02x\n", out[0], out[1], out[2], out[3], out[4], out[5], out[6], out[7]);

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
