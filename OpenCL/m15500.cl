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
#include "inc_simd.cl"
#include "inc_rp.h"
#include "inc_rp.cl"

#define COMPARE_S "inc_comp_single.cl"
#define COMPARE_M "inc_comp_multi.cl"

#define BLAKE2B_FINAL                 1
#define BLAKE2B_UPDATE                0
#define BLAKE2B_FANOUT       0x00010000
#define BLAKE2B_DEPTH        0x01000000
#define BLAKE2B_BLOCK_SIZE           64
#define ARGON2_SYNC_POINTS            4 
#define ARGON2_BLOCK_SIZE          1024
#define ARGON2_PREHASH_SEED_LENGTH   72
#define ARGON2_PREHASH_HASH_LENGTH   64

#define BLAKE2B_G(r,i,a,b,c,d)                \
  do {                                        \
    a = a + b + m[blake2b_sigma[r][2*i+0]];   \
    d = rotr64(d ^ a, 32);                    \
    c = c + d;                                \
    b = rotr64(b ^ c, 24);                    \
    a = a + b + m[blake2b_sigma[r][2*i+1]];   \
    d = rotr64(d ^ a, 16);                    \
    c = c + d;                                \
    b = rotr64(b ^ c, 63);                    \
  } while(0)

#define BLAKE2B_ROUND(r)                      \
  do {                                        \
    BLAKE2B_G(r,0,v[ 0],v[ 4],v[ 8],v[12]);   \
    BLAKE2B_G(r,1,v[ 1],v[ 5],v[ 9],v[13]);   \
    BLAKE2B_G(r,2,v[ 2],v[ 6],v[10],v[14]);   \
    BLAKE2B_G(r,3,v[ 3],v[ 7],v[11],v[15]);   \
    BLAKE2B_G(r,4,v[ 0],v[ 5],v[10],v[15]);   \
    BLAKE2B_G(r,5,v[ 1],v[ 6],v[11],v[12]);   \
    BLAKE2B_G(r,6,v[ 2],v[ 7],v[ 8],v[13]);   \
    BLAKE2B_G(r,7,v[ 3],v[ 4],v[ 9],v[14]);   \
  } while(0)

#if VECT_SIZE == 1
#define p_push( a )                              \
  do {                                           \
    p[ index++ ] = a >>  0 & 0xff;               \
    p[ index++ ] = a >>  8 & 0xff;               \
    p[ index++ ] = a >> 16 & 0xff;               \
    p[ index++ ] = a >> 24 & 0xff;               \
  } while (0)
#endif

#if VECT_SIZE == 2
#define p_push( a )                              \
  do {                                           \
    p[ index.s0++ ].s0 = a >>  0 & 0xff;         \
    p[ index.s0++ ].s0 = a >>  8 & 0xff;         \
    p[ index.s0++ ].s0 = a >> 16 & 0xff;         \
    p[ index.s0++ ].s0 = a >> 24 & 0xff;         \
    p[ index.s1++ ].s1 = a >>  0 & 0xff;         \
    p[ index.s1++ ].s1 = a >>  8 & 0xff;         \
    p[ index.s1++ ].s1 = a >> 16 & 0xff;         \
    p[ index.s1++ ].s1 = a >> 24 & 0xff;         \
  } while (0)
#endif

#if VECT_SIZE == 4
#define p_push( a )                              \
  do {                                           \
    p[ index.s0++ ].s0 = a >>  0 & 0xff;         \
    p[ index.s0++ ].s0 = a >>  8 & 0xff;         \
    p[ index.s0++ ].s0 = a >> 16 & 0xff;         \
    p[ index.s0++ ].s0 = a >> 24 & 0xff;         \
    p[ index.s1++ ].s1 = a >>  0 & 0xff;         \
    p[ index.s1++ ].s1 = a >>  8 & 0xff;         \
    p[ index.s1++ ].s1 = a >> 16 & 0xff;         \
    p[ index.s1++ ].s1 = a >> 24 & 0xff;         \
    p[ index.s2++ ].s2 = a >>  0 & 0xff;         \
    p[ index.s2++ ].s2 = a >>  8 & 0xff;         \
    p[ index.s2++ ].s2 = a >> 16 & 0xff;         \
    p[ index.s2++ ].s2 = a >> 24 & 0xff;         \
    p[ index.s3++ ].s3 = a >>  0 & 0xff;         \
    p[ index.s3++ ].s3 = a >>  8 & 0xff;         \
    p[ index.s3++ ].s3 = a >> 16 & 0xff;         \
    p[ index.s3++ ].s3 = a >> 24 & 0xff;         \
  } while (0)
#endif

#if VECT_SIZE == 8
#define p_push( a )                              \
  do {                                           \
    p[ index.s0++ ].s0 = a >>  0 & 0xff;         \
    p[ index.s0++ ].s0 = a >>  8 & 0xff;         \
    p[ index.s0++ ].s0 = a >> 16 & 0xff;         \
    p[ index.s0++ ].s0 = a >> 24 & 0xff;         \
    p[ index.s1++ ].s1 = a >>  0 & 0xff;         \
    p[ index.s1++ ].s1 = a >>  8 & 0xff;         \
    p[ index.s1++ ].s1 = a >> 16 & 0xff;         \
    p[ index.s1++ ].s1 = a >> 24 & 0xff;         \
    p[ index.s2++ ].s2 = a >>  0 & 0xff;         \
    p[ index.s2++ ].s2 = a >>  8 & 0xff;         \
    p[ index.s2++ ].s2 = a >> 16 & 0xff;         \
    p[ index.s2++ ].s2 = a >> 24 & 0xff;         \
    p[ index.s3++ ].s3 = a >>  0 & 0xff;         \
    p[ index.s3++ ].s3 = a >>  8 & 0xff;         \
    p[ index.s3++ ].s3 = a >> 16 & 0xff;         \
    p[ index.s3++ ].s3 = a >> 24 & 0xff;         \
    p[ index.s4++ ].s4 = a >>  0 & 0xff;         \
    p[ index.s4++ ].s4 = a >>  8 & 0xff;         \
    p[ index.s4++ ].s4 = a >> 16 & 0xff;         \
    p[ index.s4++ ].s4 = a >> 24 & 0xff;         \
    p[ index.s5++ ].s5 = a >>  0 & 0xff;         \
    p[ index.s5++ ].s5 = a >>  8 & 0xff;         \
    p[ index.s5++ ].s5 = a >> 16 & 0xff;         \
    p[ index.s5++ ].s5 = a >> 24 & 0xff;         \
    p[ index.s6++ ].s6 = a >>  0 & 0xff;         \
    p[ index.s6++ ].s6 = a >>  8 & 0xff;         \
    p[ index.s6++ ].s6 = a >> 16 & 0xff;         \
    p[ index.s6++ ].s6 = a >> 24 & 0xff;         \
    p[ index.s7++ ].s7 = a >>  0 & 0xff;         \
    p[ index.s7++ ].s7 = a >>  8 & 0xff;         \
    p[ index.s7++ ].s7 = a >> 16 & 0xff;         \
    p[ index.s7++ ].s7 = a >> 24 & 0xff;         \
  } while (0)
#endif

#if VECT_SIZE == 16
#define p_push( a )                              \
  do {                                           \
    p[ index.s0++ ].s0 = a >>  0 & 0xff;         \
    p[ index.s0++ ].s0 = a >>  8 & 0xff;         \
    p[ index.s0++ ].s0 = a >> 16 & 0xff;         \
    p[ index.s0++ ].s0 = a >> 24 & 0xff;         \
    p[ index.s1++ ].s1 = a >>  0 & 0xff;         \
    p[ index.s1++ ].s1 = a >>  8 & 0xff;         \
    p[ index.s1++ ].s1 = a >> 16 & 0xff;         \
    p[ index.s1++ ].s1 = a >> 24 & 0xff;         \
    p[ index.s2++ ].s2 = a >>  0 & 0xff;         \
    p[ index.s2++ ].s2 = a >>  8 & 0xff;         \
    p[ index.s2++ ].s2 = a >> 16 & 0xff;         \
    p[ index.s2++ ].s2 = a >> 24 & 0xff;         \
    p[ index.s3++ ].s3 = a >>  0 & 0xff;         \
    p[ index.s3++ ].s3 = a >>  8 & 0xff;         \
    p[ index.s3++ ].s3 = a >> 16 & 0xff;         \
    p[ index.s3++ ].s3 = a >> 24 & 0xff;         \
    p[ index.s4++ ].s4 = a >>  0 & 0xff;         \
    p[ index.s4++ ].s4 = a >>  8 & 0xff;         \
    p[ index.s4++ ].s4 = a >> 16 & 0xff;         \
    p[ index.s4++ ].s4 = a >> 24 & 0xff;         \
    p[ index.s5++ ].s5 = a >>  0 & 0xff;         \
    p[ index.s5++ ].s5 = a >>  8 & 0xff;         \
    p[ index.s5++ ].s5 = a >> 16 & 0xff;         \
    p[ index.s5++ ].s5 = a >> 24 & 0xff;         \
    p[ index.s6++ ].s6 = a >>  0 & 0xff;         \
    p[ index.s6++ ].s6 = a >>  8 & 0xff;         \
    p[ index.s6++ ].s6 = a >> 16 & 0xff;         \
    p[ index.s6++ ].s6 = a >> 24 & 0xff;         \
    p[ index.s7++ ].s7 = a >>  0 & 0xff;         \
    p[ index.s7++ ].s7 = a >>  8 & 0xff;         \
    p[ index.s7++ ].s7 = a >> 16 & 0xff;         \
    p[ index.s7++ ].s7 = a >> 24 & 0xff;         \
    p[ index.s8++ ].s8 = a >>  0 & 0xff;         \
    p[ index.s8++ ].s8 = a >>  8 & 0xff;         \
    p[ index.s8++ ].s8 = a >> 16 & 0xff;         \
    p[ index.s8++ ].s8 = a >> 24 & 0xff;         \
    p[ index.s9++ ].s9 = a >>  0 & 0xff;         \
    p[ index.s9++ ].s9 = a >>  8 & 0xff;         \
    p[ index.s9++ ].s9 = a >> 16 & 0xff;         \
    p[ index.s9++ ].s9 = a >> 24 & 0xff;         \
    p[ index.sa++ ].sa = a >>  0 & 0xff;         \
    p[ index.sa++ ].sa = a >>  8 & 0xff;         \
    p[ index.sa++ ].sa = a >> 16 & 0xff;         \
    p[ index.sa++ ].sa = a >> 24 & 0xff;         \
    p[ index.sb++ ].sb = a >>  0 & 0xff;         \
    p[ index.sb++ ].sb = a >>  8 & 0xff;         \
    p[ index.sb++ ].sb = a >> 16 & 0xff;         \
    p[ index.sb++ ].sb = a >> 24 & 0xff;         \
    p[ index.sc++ ].sc = a >>  0 & 0xff;         \
    p[ index.sc++ ].sc = a >>  8 & 0xff;         \
    p[ index.sc++ ].sc = a >> 16 & 0xff;         \
    p[ index.sc++ ].sc = a >> 24 & 0xff;         \
    p[ index.sd++ ].sd = a >>  0 & 0xff;         \
    p[ index.sd++ ].sd = a >>  8 & 0xff;         \
    p[ index.sd++ ].sd = a >> 16 & 0xff;         \
    p[ index.sd++ ].sd = a >> 24 & 0xff;         \
    p[ index.se++ ].se = a >>  0 & 0xff;         \
    p[ index.se++ ].se = a >>  8 & 0xff;         \
    p[ index.se++ ].se = a >> 16 & 0xff;         \
    p[ index.se++ ].se = a >> 24 & 0xff;         \
    p[ index.sf++ ].sf = a >>  0 & 0xff;         \
    p[ index.sf++ ].sf = a >>  8 & 0xff;         \
    p[ index.sf++ ].sf = a >> 16 & 0xff;         \
    p[ index.sf++ ].sf = a >> 24 & 0xff;         \
  } while (0)
#endif

inline u32 get_max(const u32x a)
{
  u32 r = 0;

  #if VECT_SIZE == 1
  r = a;
  #endif

  #if VECT_SIZE >= 2
  r = r < a.s0 ? a.s0 : r;
  r = r < a.s1 ? a.s1 : r;
  #endif

  #if VECT_SIZE >= 4
  r = r < a.s2 ? a.s2 : r;
  r = r < a.s3 ? a.s3 : r;
  #endif

  #if VECT_SIZE >= 8
  r = r < a.s4 ? a.s4 : r;
  r = r < a.s5 ? a.s5 : r;
  r = r < a.s6 ? a.s6 : r;
  r = r < a.s7 ? a.s7 : r;
  #endif

  #if VECT_SIZE >=16
  r = r < a.s8 ? a.s8 : r;
  r = r < a.s9 ? a.s9 : r;
  r = r < a.sa ? a.sa : r;
  r = r < a.sb ? a.sb : r;
  r = r < a.sc ? a.sc : r;
  r = r < a.sd ? a.sd : r;
  r = r < a.se ? a.se : r;
  r = r < a.sf ? a.sf : r;
  #endif

  return r;
}

void blake2b_transform(u64x h[8], u64x t[2], u64x f[2], u64x m[16], u64x v[16], const u64x out_len, const u8 isFinal)
{
  if (isFinal)
    f[0] = -1;

  t[0] += out_len;

  v[ 0] = h[0];
  v[ 1] = h[1];
  v[ 2] = h[2];
  v[ 3] = h[3];
  v[ 4] = h[4];
  v[ 5] = h[5];
  v[ 6] = h[6];
  v[ 7] = h[7];
  v[ 8] = BLAKE2B_IV_00;
  v[ 9] = BLAKE2B_IV_01;
  v[10] = BLAKE2B_IV_02;
  v[11] = BLAKE2B_IV_03;
  v[12] = BLAKE2B_IV_04 ^ t[0];
  v[13] = BLAKE2B_IV_05 ^ t[1];
  v[14] = BLAKE2B_IV_06 ^ f[0];
  v[15] = BLAKE2B_IV_07 ^ f[1];

  const u8a blake2b_sigma[12][16] =
  {
    {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15 } ,
    { 14, 10,  4,  8,  9, 15, 13,  6,  1, 12,  0,  2, 11,  7,  5,  3 } ,
    { 11,  8, 12,  0,  5,  2, 15, 13, 10, 14,  3,  6,  7,  1,  9,  4 } ,
    {  7,  9,  3,  1, 13, 12, 11, 14,  2,  6,  5, 10,  4,  0, 15,  8 } ,
    {  9,  0,  5,  7,  2,  4, 10, 15, 14,  1, 11, 12,  6,  8,  3, 13 } ,
    {  2, 12,  6, 10,  0, 11,  8,  3,  4, 13,  7,  5, 15, 14,  1,  9 } ,
    { 12,  5,  1, 15, 14, 13,  4, 10,  0,  7,  6,  3,  9,  2,  8, 11 } ,
    { 13, 11,  7, 14, 12,  1,  3,  9,  5,  0, 15,  4,  8,  6,  2, 10 } ,
    {  6, 15, 14,  9, 11,  3,  0,  8, 12,  2, 13,  7,  1,  4, 10,  5 } ,
    { 10,  2,  8,  4,  7,  6,  1,  5, 15, 11,  9, 14,  3, 12, 13 , 0 } ,
    {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15 } ,
    { 14, 10,  4,  8,  9, 15, 13,  6,  1, 12,  0,  2, 11,  7,  5,  3 }
  };

  BLAKE2B_ROUND( 0);
  BLAKE2B_ROUND( 1);
  BLAKE2B_ROUND( 2);
  BLAKE2B_ROUND( 3);
  BLAKE2B_ROUND( 4);
  BLAKE2B_ROUND( 5);
  BLAKE2B_ROUND( 6);
  BLAKE2B_ROUND( 7);
  BLAKE2B_ROUND( 8);
  BLAKE2B_ROUND( 9);
  BLAKE2B_ROUND(10);
  BLAKE2B_ROUND(11);

  h[0] = h[0] ^ v[0] ^ v[ 8];
  h[1] = h[1] ^ v[1] ^ v[ 9];
  h[2] = h[2] ^ v[2] ^ v[10];
  h[3] = h[3] ^ v[3] ^ v[11];
  h[4] = h[4] ^ v[4] ^ v[12];
  h[5] = h[5] ^ v[5] ^ v[13];
  h[6] = h[6] ^ v[6] ^ v[14];
  h[7] = h[7] ^ v[7] ^ v[15];

}

__kernel void m15500_init (__global pw_t *pws, __global const kernel_rule_t *rules_buf, __global const comb_t *combs_buf, __global const bf_t *bfs_buf, __global argon2_tmp_t *tmps, __global void *hooks, __global const u32 *bitmaps_buf_s1_a, __global const u32 *bitmaps_buf_s1_b, __global const u32 *bitmaps_buf_s1_c, __global const u32 *bitmaps_buf_s1_d, __global const u32 *bitmaps_buf_s2_a, __global const u32 *bitmaps_buf_s2_b, __global const u32 *bitmaps_buf_s2_c, __global const u32 *bitmaps_buf_s2_d, __global plain_t *plains_buf, __global const digest_t *digests_buf, __global u32 *hashes_shown, __global const salt_t *salt_bufs, __global argon2_t *esalt_bufs, __global u32 *d_return_buf, __global u32 *d_scryptV0_buf, __global u32 *d_scryptV1_buf, __global u32 *d_scryptV2_buf, __global u32 *d_scryptV3_buf, const u32 bitmap_mask, const u32 bitmap_shift1, const u32 bitmap_shift2, const u32 salt_pos, const u32 loop_pos, const u32 loop_cnt, const u32 rules_cnt, const u32 digests_cnt, const u32 digests_offset, const u32 combs_mode, const u32 gid_max)
{
  /**
   * base
   */

  const u32 gid   = get_global_id (0);

  salt_t   salt   = salt_bufs[salt_pos];
  argon2_t argon2 = esalt_bufs[salt_pos];

  if (gid >= gid_max) return;

  u32x w0[8];

  w0[0] = pws[gid].i[0];
  w0[1] = pws[gid].i[1];
  w0[2] = pws[gid].i[2];
  w0[3] = pws[gid].i[3];
  w0[4] = pws[gid].i[4];
  w0[5] = pws[gid].i[5];
  w0[6] = pws[gid].i[6];
  w0[7] = pws[gid].i[7];

  const u32x pw_len = pws[gid].pw_len;

  /**
   * Fill byte[] buffer with Argon2 Parameters 
   * 
   * H (p, l, m, t, v, y, <P>, P, <S>, S, <K>, K, <X>, X)
   */

  u64x   m[ 16] = { 0 };
  u64x   v[ 16] = { 0 };
  u64x   h[  8] = { 0 };
  u64x   f[  2] = { 0 };
  u64x   t[  2] = { 0 };
  u64x pre[  8] = { 0 };
  u8x    p[128] = { 0 };
  u32x   index  =   0  ; 

  p_push(argon2.p);                                     // p
  p_push(argon2.l);                                     // l
  p_push(argon2.m);                                     // m
  p_push(argon2.t);                                     // t
  p_push(argon2.v);                                     // v
  p_push(argon2.y);                                     // y

                                                        // <P>
  #if VECT_SIZE == 1
    p[ index++ ] = pw_len >>  0 & 0xff;
    p[ index++ ] = pw_len >>  8 & 0xff;
    p[ index++ ] = pw_len >> 16 & 0xff;
    p[ index++ ] = pw_len >> 24 & 0xff;
  #endif

  #if VECT_SIZE >= 2
    p[ index.s0++ ].s0 = pw_len.s0 >>  0 & 0xff;
    p[ index.s0++ ].s0 = pw_len.s0 >>  8 & 0xff;
    p[ index.s0++ ].s0 = pw_len.s0 >> 16 & 0xff;
    p[ index.s0++ ].s0 = pw_len.s0 >> 24 & 0xff;
    p[ index.s1++ ].s1 = pw_len.s1 >>  0 & 0xff;
    p[ index.s1++ ].s1 = pw_len.s1 >>  8 & 0xff;
    p[ index.s1++ ].s1 = pw_len.s1 >> 16 & 0xff;
    p[ index.s1++ ].s1 = pw_len.s1 >> 24 & 0xff;
  #endif

  #if VECT_SIZE >= 4
    p[ index.s2++ ].s2 = pw_len.s2 >>  0 & 0xff;
    p[ index.s2++ ].s2 = pw_len.s2 >>  8 & 0xff;
    p[ index.s2++ ].s2 = pw_len.s2 >> 16 & 0xff;
    p[ index.s2++ ].s2 = pw_len.s2 >> 24 & 0xff;
    p[ index.s3++ ].s3 = pw_len.s3 >>  0 & 0xff;
    p[ index.s3++ ].s3 = pw_len.s3 >>  8 & 0xff;
    p[ index.s3++ ].s3 = pw_len.s3 >> 16 & 0xff;
    p[ index.s3++ ].s3 = pw_len.s3 >> 24 & 0xff;
  #endif
  
  #if VECT_SIZE >= 8
    p[ index.s4++ ].s4 = pw_len.s4 >>  0 & 0xff;
    p[ index.s4++ ].s4 = pw_len.s4 >>  8 & 0xff;
    p[ index.s4++ ].s4 = pw_len.s4 >> 16 & 0xff;
    p[ index.s4++ ].s4 = pw_len.s4 >> 24 & 0xff;
    p[ index.s5++ ].s5 = pw_len.s5 >>  0 & 0xff;
    p[ index.s5++ ].s5 = pw_len.s5 >>  8 & 0xff;
    p[ index.s5++ ].s5 = pw_len.s5 >> 16 & 0xff;
    p[ index.s5++ ].s5 = pw_len.s5 >> 24 & 0xff;
    p[ index.s6++ ].s6 = pw_len.s6 >>  0 & 0xff;
    p[ index.s6++ ].s6 = pw_len.s6 >>  8 & 0xff;
    p[ index.s6++ ].s6 = pw_len.s6 >> 16 & 0xff;
    p[ index.s6++ ].s6 = pw_len.s6 >> 24 & 0xff;
    p[ index.s7++ ].s7 = pw_len.s7 >>  0 & 0xff;
    p[ index.s7++ ].s7 = pw_len.s7 >>  8 & 0xff;
    p[ index.s7++ ].s7 = pw_len.s7 >> 16 & 0xff;
    p[ index.s7++ ].s7 = pw_len.s7 >> 24 & 0xff;
  #endif

  #if VECT_SIZE >= 16
    p[ index.s8++ ].s8 = pw_len.s8 >>  0 & 0xff;
    p[ index.s8++ ].s8 = pw_len.s8 >>  8 & 0xff;
    p[ index.s8++ ].s8 = pw_len.s8 >> 16 & 0xff;
    p[ index.s8++ ].s8 = pw_len.s8 >> 24 & 0xff;
    p[ index.s9++ ].s9 = pw_len.s9 >>  0 & 0xff;
    p[ index.s9++ ].s9 = pw_len.s9 >>  8 & 0xff;
    p[ index.s9++ ].s9 = pw_len.s9 >> 16 & 0xff;
    p[ index.s9++ ].s9 = pw_len.s9 >> 24 & 0xff;
    p[ index.sa++ ].sa = pw_len.sa >>  0 & 0xff;
    p[ index.sa++ ].sa = pw_len.sa >>  8 & 0xff;
    p[ index.sa++ ].sa = pw_len.sa >> 16 & 0xff;
    p[ index.sa++ ].sa = pw_len.sa >> 24 & 0xff;
    p[ index.sb++ ].sb = pw_len.sb >>  0 & 0xff;
    p[ index.sb++ ].sb = pw_len.sb >>  8 & 0xff;
    p[ index.sb++ ].sb = pw_len.sb >> 16 & 0xff;
    p[ index.sb++ ].sb = pw_len.sb >> 24 & 0xff;
    p[ index.sc++ ].sc = pw_len.sc >>  0 & 0xff;
    p[ index.sc++ ].sc = pw_len.sc >>  8 & 0xff;
    p[ index.sc++ ].sc = pw_len.sc >> 16 & 0xff;
    p[ index.sc++ ].sc = pw_len.sc >> 24 & 0xff;
    p[ index.sd++ ].sd = pw_len.sd >>  0 & 0xff;
    p[ index.sd++ ].sd = pw_len.sd >>  8 & 0xff;
    p[ index.sd++ ].sd = pw_len.sd >> 16 & 0xff;
    p[ index.sd++ ].sd = pw_len.sd >> 24 & 0xff;
    p[ index.se++ ].se = pw_len.se >>  0 & 0xff;
    p[ index.se++ ].se = pw_len.se >>  8 & 0xff;
    p[ index.se++ ].se = pw_len.se >> 16 & 0xff;
    p[ index.se++ ].se = pw_len.se >> 24 & 0xff;
    p[ index.sf++ ].sf = pw_len.sf >>  0 & 0xff;
    p[ index.sf++ ].sf = pw_len.sf >>  8 & 0xff;
    p[ index.sf++ ].sf = pw_len.sf >> 16 & 0xff;
    p[ index.sf++ ].sf = pw_len.sf >> 24 & 0xff;
  #endif

                                                                                   // P
  for (u8 i = 0; i < get_max(pw_len); i++)
  {
    #if VECT_SIZE == 1
      p[ index++ ] = convert_uchar(w0[i / 4] >> ((i % 4) * 8) & 0xff);
    #endif 

    #if VECT_SIZE >= 2
      p[ index.s0++ ].s0 = convert_uchar(w0[i / 4].s0 >> ((i % 4) * 8) & 0xff);
      p[ index.s1++ ].s1 = convert_uchar(w0[i / 4].s1 >> ((i % 4) * 8) & 0xff);
    #endif

    #if VECT_SIZE >= 4
      p[ index.s2++ ].s2 = convert_uchar(w0[i / 4].s2 >> ((i % 4) * 8) & 0xff);
      p[ index.s3++ ].s3 = convert_uchar(w0[i / 4].s3 >> ((i % 4) * 8) & 0xff);
    #endif

    #if VECT_SIZE >= 8
      p[ index.s4++ ].s4 = convert_uchar(w0[i / 4].s4 >> ((i % 4) * 8) & 0xff);
      p[ index.s5++ ].s5 = convert_uchar(w0[i / 4].s5 >> ((i % 4) * 8) & 0xff);
      p[ index.s6++ ].s6 = convert_uchar(w0[i / 4].s6 >> ((i % 4) * 8) & 0xff);
      p[ index.s7++ ].s7 = convert_uchar(w0[i / 4].s7 >> ((i % 4) * 8) & 0xff);
    #endif

    #if VECT_SIZE >= 16
      p[ index.s8++ ].s8 = convert_uchar(w0[i / 4].s8 >> ((i % 4) * 8) & 0xff);
      p[ index.s9++ ].s9 = convert_uchar(w0[i / 4].s9 >> ((i % 4) * 8) & 0xff);
      p[ index.sa++ ].sa = convert_uchar(w0[i / 4].sa >> ((i % 4) * 8) & 0xff);
      p[ index.sb++ ].sb = convert_uchar(w0[i / 4].sb >> ((i % 4) * 8) & 0xff);
      p[ index.sc++ ].sc = convert_uchar(w0[i / 4].sc >> ((i % 4) * 8) & 0xff);
      p[ index.sd++ ].sd = convert_uchar(w0[i / 4].sd >> ((i % 4) * 8) & 0xff);
      p[ index.se++ ].se = convert_uchar(w0[i / 4].se >> ((i % 4) * 8) & 0xff);
      p[ index.sf++ ].sf = convert_uchar(w0[i / 4].sf >> ((i % 4) * 8) & 0xff);
    #endif
  }

  p_push(salt.salt_len);                                                           // <S>

  for (u8 i = 0; i < get_max(salt.salt_len); i++)                                  //  S
  { 
    #if VECT_SIZE == 1 
      p[ index++ ] = salt.salt_buf[i / 4] >> ((i % 4) * 8) & 0xff;
    #endif

    #if VECT_SIZE >= 2
      p[ index.s0++ ].s0 = salt.salt_buf[i / 4] >> ((i % 4) * 8) & 0xff; 
      p[ index.s1++ ].s1 = salt.salt_buf[i / 4] >> ((i % 4) * 8) & 0xff;
    #endif

    #if VECT_SIZE >= 4
      p[ index.s2++ ].s2 = salt.salt_buf[i / 4] >> ((i % 4) * 8) & 0xff;
      p[ index.s3++ ].s3 = salt.salt_buf[i / 4] >> ((i % 4) * 8) & 0xff;
    #endif

    #if VECT_SIZE >= 8
      p[ index.s4++ ].s4 = salt.salt_buf[i / 4] >> ((i % 4) * 8) & 0xff;
      p[ index.s5++ ].s5 = salt.salt_buf[i / 4] >> ((i % 4) * 8) & 0xff;
      p[ index.s6++ ].s6 = salt.salt_buf[i / 4] >> ((i % 4) * 8) & 0xff;
      p[ index.s7++ ].s7 = salt.salt_buf[i / 4] >> ((i % 4) * 8) & 0xff;
    #endif

    #if VECT_SIZE >= 16
      p[ index.s8++ ].s8 = salt.salt_buf[i / 4] >> ((i % 4) * 8) & 0xff;
      p[ index.s9++ ].s9 = salt.salt_buf[i / 4] >> ((i % 4) * 8) & 0xff;
      p[ index.sa++ ].sa = salt.salt_buf[i / 4] >> ((i % 4) * 8) & 0xff;
      p[ index.sb++ ].sb = salt.salt_buf[i / 4] >> ((i % 4) * 8) & 0xff;
      p[ index.sc++ ].sc = salt.salt_buf[i / 4] >> ((i % 4) * 8) & 0xff;
      p[ index.sd++ ].sd = salt.salt_buf[i / 4] >> ((i % 4) * 8) & 0xff;
      p[ index.se++ ].se = salt.salt_buf[i / 4] >> ((i % 4) * 8) & 0xff;
      p[ index.sf++ ].sf = salt.salt_buf[i / 4] >> ((i % 4) * 8) & 0xff;
    #endif
  }

  /*
   * Convert byte[] buffer to 64 bit vectors for Blake2b compress function
   */

  for (u8 i = 0; i < 16; i++)
  {
    #if VECT_SIZE == 1
      m[i]  ^=  p[i * 8 + 4] <<  0;
      m[i]  ^=  p[i * 8 + 5] <<  8;
      m[i]  ^=  p[i * 8 + 6] << 16;
      m[i]  ^=  p[i * 8 + 7] << 24;
      m[i] <<=  32;
      m[i]  ^=  p[i * 8 + 0] <<  0;
      m[i]  ^=  p[i * 8 + 1] <<  8;
      m[i]  ^=  p[i * 8 + 2] << 16;
      m[i]  ^=  p[i * 8 + 3] << 24;
    #endif  

    #if VECT_SIZE >= 2
      m[i].s0  ^=  p[i * 8 + 4].s0 <<  0;
      m[i].s0  ^=  p[i * 8 + 5].s0 <<  8;
      m[i].s0  ^=  p[i * 8 + 6].s0 << 16;
      m[i].s0  ^=  p[i * 8 + 7].s0 << 24;
      m[i].s0 <<=                     32;
      m[i].s0  ^=  p[i * 8 + 0].s0 <<  0;
      m[i].s0  ^=  p[i * 8 + 1].s0 <<  8;
      m[i].s0  ^=  p[i * 8 + 2].s0 << 16;
      m[i].s0  ^=  p[i * 8 + 3].s0 << 24;

      m[i].s1  ^=  p[i * 8 + 4].s1 <<  0;
      m[i].s1  ^=  p[i * 8 + 5].s1 <<  8;
      m[i].s1  ^=  p[i * 8 + 6].s1 << 16;
      m[i].s1  ^=  p[i * 8 + 7].s1 << 24;
      m[i].s1 <<=                     32;
      m[i].s1  ^=  p[i * 8 + 0].s1 <<  0;
      m[i].s1  ^=  p[i * 8 + 1].s1 <<  8;
      m[i].s1  ^=  p[i * 8 + 2].s1 << 16;
      m[i].s1  ^=  p[i * 8 + 3].s1 << 24;
    #endif

    #if VECT_SIZE >= 4
      m[i].s2  ^=  p[i * 8 + 4].s2 <<  0;
      m[i].s2  ^=  p[i * 8 + 5].s2 <<  8;
      m[i].s2  ^=  p[i * 8 + 6].s2 << 16;
      m[i].s2  ^=  p[i * 8 + 7].s2 << 24;
      m[i].s2 <<=                     32;
      m[i].s2  ^=  p[i * 8 + 0].s2 <<  0;
      m[i].s2  ^=  p[i * 8 + 1].s2 <<  8;
      m[i].s2  ^=  p[i * 8 + 2].s2 << 16;
      m[i].s2  ^=  p[i * 8 + 3].s2 << 24;

      m[i].s3  ^=  p[i * 8 + 4].s2 <<  0;
      m[i].s3  ^=  p[i * 8 + 5].s2 <<  8;
      m[i].s3  ^=  p[i * 8 + 6].s2 << 16;
      m[i].s3  ^=  p[i * 8 + 7].s2 << 24;
      m[i].s3 <<=                     32;
      m[i].s3  ^=  p[i * 8 + 0].s2 <<  0;
      m[i].s3  ^=  p[i * 8 + 1].s2 <<  8;
      m[i].s3  ^=  p[i * 8 + 2].s2 << 16;
      m[i].s3  ^=  p[i * 8 + 3].s2 << 24;
    #endif

    #if VECT_SIZE >= 8
      m[i].s4  ^=  p[i * 8 + 4].s4 <<  0;
      m[i].s4  ^=  p[i * 8 + 5].s4 <<  8;
      m[i].s4  ^=  p[i * 8 + 6].s4 << 16;
      m[i].s4  ^=  p[i * 8 + 7].s4 << 24;
      m[i].s4 <<=                     32;
      m[i].s4  ^=  p[i * 8 + 0].s4 <<  0;
      m[i].s4  ^=  p[i * 8 + 1].s4 <<  8;
      m[i].s4  ^=  p[i * 8 + 2].s4 << 16;
      m[i].s4  ^=  p[i * 8 + 3].s4 << 24;

      m[i].s5  ^=  p[i * 8 + 4].s5 <<  0;
      m[i].s5  ^=  p[i * 8 + 5].s5 <<  8;
      m[i].s5  ^=  p[i * 8 + 6].s5 << 16;
      m[i].s5  ^=  p[i * 8 + 7].s5 << 24;
      m[i].s5 <<=                     32;
      m[i].s5  ^=  p[i * 8 + 0].s5 <<  0;
      m[i].s5  ^=  p[i * 8 + 1].s5 <<  8;
      m[i].s5  ^=  p[i * 8 + 2].s5 << 16;
      m[i].s5  ^=  p[i * 8 + 3].s5 << 24;

      m[i].s6  ^=  p[i * 8 + 4].s6 <<  0;
      m[i].s6  ^=  p[i * 8 + 5].s6 <<  8;
      m[i].s6  ^=  p[i * 8 + 6].s6 << 16;
      m[i].s6  ^=  p[i * 8 + 7].s6 << 24;
      m[i].s6 <<=                     32;
      m[i].s6  ^=  p[i * 8 + 0].s6 <<  0;
      m[i].s6  ^=  p[i * 8 + 1].s6 <<  8;
      m[i].s6  ^=  p[i * 8 + 2].s6 << 16;
      m[i].s6  ^=  p[i * 8 + 3].s6 << 24;

      m[i].s7  ^=  p[i * 8 + 4].s7 <<  0;
      m[i].s7  ^=  p[i * 8 + 5].s7 <<  8;
      m[i].s7  ^=  p[i * 8 + 6].s7 << 16;
      m[i].s7  ^=  p[i * 8 + 7].s7 << 24;
      m[i].s7 <<=                     32;
      m[i].s7  ^=  p[i * 8 + 0].s7 <<  0;
      m[i].s7  ^=  p[i * 8 + 1].s7 <<  8;
      m[i].s7  ^=  p[i * 8 + 2].s7 << 16;
      m[i].s7  ^=  p[i * 8 + 3].s7 << 24;
    #endif

    #if VECT_SIZE >= 16
      m[i].s8  ^=  p[i * 8 + 4].s8 <<  0;
      m[i].s8  ^=  p[i * 8 + 5].s8 <<  8;
      m[i].s8  ^=  p[i * 8 + 6].s8 << 16;
      m[i].s8  ^=  p[i * 8 + 7].s8 << 24;
      m[i].s8 <<=                     32;
      m[i].s8  ^=  p[i * 8 + 0].s8 <<  0;
      m[i].s8  ^=  p[i * 8 + 1].s8 <<  8;
      m[i].s8  ^=  p[i * 8 + 2].s8 << 16;
      m[i].s8  ^=  p[i * 8 + 3].s8 << 24;

      m[i].s9  ^=  p[i * 8 + 4].s9 <<  0;
      m[i].s9  ^=  p[i * 8 + 5].s9 <<  8;
      m[i].s9  ^=  p[i * 8 + 6].s9 << 16;
      m[i].s9  ^=  p[i * 8 + 7].s9 << 24;
      m[i].s9 <<=                     32;
      m[i].s9  ^=  p[i * 8 + 0].s9 <<  0;
      m[i].s9  ^=  p[i * 8 + 1].s9 <<  8;
      m[i].s9  ^=  p[i * 8 + 2].s9 << 16;
      m[i].s9  ^=  p[i * 8 + 3].s9 << 24;

      m[i].sa  ^=  p[i * 8 + 4].sa <<  0;
      m[i].sa  ^=  p[i * 8 + 5].sa <<  8;
      m[i].sa  ^=  p[i * 8 + 6].sa << 16;
      m[i].sa  ^=  p[i * 8 + 7].sa << 24;
      m[i].sa <<=                     32;
      m[i].sa  ^=  p[i * 8 + 0].sa <<  0;
      m[i].sa  ^=  p[i * 8 + 1].sa <<  8;
      m[i].sa  ^=  p[i * 8 + 2].sa << 16;
      m[i].sa  ^=  p[i * 8 + 3].sa << 24;

      m[i].sb  ^=  p[i * 8 + 4].sb <<  0;
      m[i].sb  ^=  p[i * 8 + 5].sb <<  8;
      m[i].sb  ^=  p[i * 8 + 6].sb << 16;
      m[i].sb  ^=  p[i * 8 + 7].sb << 24;
      m[i].sb <<=                     32;
      m[i].sb  ^=  p[i * 8 + 0].sb <<  0;
      m[i].sb  ^=  p[i * 8 + 1].sb <<  8;
      m[i].sb  ^=  p[i * 8 + 2].sb << 16;
      m[i].sb  ^=  p[i * 8 + 3].sb << 24;

      m[i].sc  ^=  p[i * 8 + 4].sc <<  0;
      m[i].sc  ^=  p[i * 8 + 5].sc <<  8;
      m[i].sc  ^=  p[i * 8 + 6].sc << 16;
      m[i].sc  ^=  p[i * 8 + 7].sc << 24;
      m[i].sc <<=                     32;
      m[i].sc  ^=  p[i * 8 + 0].sc <<  0;
      m[i].sc  ^=  p[i * 8 + 1].sc <<  8;
      m[i].sc  ^=  p[i * 8 + 2].sc << 16;
      m[i].sc  ^=  p[i * 8 + 3].sc << 24;

      m[i].sd  ^=  p[i * 8 + 4].sd <<  0;
      m[i].sd  ^=  p[i * 8 + 5].sd <<  8;
      m[i].sd  ^=  p[i * 8 + 6].sd << 16;
      m[i].sd  ^=  p[i * 8 + 7].sd << 24;
      m[i].sd <<=                     32;
      m[i].sd  ^=  p[i * 8 + 0].sd <<  0;
      m[i].sd  ^=  p[i * 8 + 1].sd <<  8;
      m[i].sd  ^=  p[i * 8 + 2].sd << 16;
      m[i].sd  ^=  p[i * 8 + 3].sd << 24;

      m[i].se  ^=  p[i * 8 + 4].se <<  0;
      m[i].se  ^=  p[i * 8 + 5].se <<  8;
      m[i].se  ^=  p[i * 8 + 6].se << 16;
      m[i].se  ^=  p[i * 8 + 7].se << 24;
      m[i].se <<=                     32;
      m[i].se  ^=  p[i * 8 + 0].se <<  0;
      m[i].se  ^=  p[i * 8 + 1].se <<  8;
      m[i].se  ^=  p[i * 8 + 2].se << 16;
      m[i].se  ^=  p[i * 8 + 3].se << 24;

      m[i].sf  ^=  p[i * 8 + 4].sf <<  0;
      m[i].sf  ^=  p[i * 8 + 5].sf <<  8;
      m[i].sf  ^=  p[i * 8 + 6].sf << 16;
      m[i].sf  ^=  p[i * 8 + 7].sf << 24;
      m[i].sf <<=                     32;
      m[i].sf  ^=  p[i * 8 + 0].sf <<  0;
      m[i].sf  ^=  p[i * 8 + 1].sf <<  8;
      m[i].sf  ^=  p[i * 8 + 2].sf << 16;
      m[i].sf  ^=  p[i * 8 + 3].sf << 24;

    #endif
  }
  
  // <K> is set to zero
  //  K  is dropped
  // <X> is set to zero
  //  X  is dropped

  /*
   * Init0 Blake2b
   */

  h[0]  = BLAKE2B_IV_00; 
  h[0] ^= BLAKE2B_BLOCK_SIZE;
  h[0] ^= BLAKE2B_FANOUT;
  h[0] ^= BLAKE2B_DEPTH;
  h[1]  = BLAKE2B_IV_01;
  h[2]  = BLAKE2B_IV_02;
  h[3]  = BLAKE2B_IV_03;
  h[4]  = BLAKE2B_IV_04;
  h[5]  = BLAKE2B_IV_05;
  h[6]  = BLAKE2B_IV_06;
  h[7]  = BLAKE2B_IV_07;

  u64x size_in = hl32_to_64(0, pw_len) + salt.salt_len + 40;

  /*
   * Blake2b transform to obtain the pre-computed hash
   */

  blake2b_transform(h, t, f, m, v, size_in, BLAKE2B_FINAL);

  pre[0] = h[0];
  pre[1] = h[1];
  pre[2] = h[2];
  pre[3] = h[3];
  pre[4] = h[4];
  pre[5] = h[5];
  pre[6] = h[6];
  pre[7] = h[7];

  /*
   * Compute two memory blocks per lane with pre-computed hash
   */

  for (u32 lane = 0; lane < argon2.p; lane++)
  {
    h[0]  = BLAKE2B_IV_00;
    h[0] ^= BLAKE2B_BLOCK_SIZE;
    h[0] ^= BLAKE2B_FANOUT;
    h[0] ^= BLAKE2B_DEPTH;
    h[1]  = BLAKE2B_IV_01;
    h[2]  = BLAKE2B_IV_02;
    h[3]  = BLAKE2B_IV_03;
    h[4]  = BLAKE2B_IV_04;
    h[5]  = BLAKE2B_IV_05;
    h[6]  = BLAKE2B_IV_06;
    h[7]  = BLAKE2B_IV_07;  

    f[0] = 0;
    f[1] = 0;
    t[0] = 0;
    t[1] = 0;

    m[0] = hl32_to_64(l32_from_64(pre[0]), ARGON2_BLOCK_SIZE);
    m[1] = hl32_to_64(l32_from_64(pre[1]), h32_from_64(pre[0]));
    m[2] = hl32_to_64(l32_from_64(pre[2]), h32_from_64(pre[1]));
    m[3] = hl32_to_64(l32_from_64(pre[3]), h32_from_64(pre[2]));
    m[4] = hl32_to_64(l32_from_64(pre[4]), h32_from_64(pre[3]));
    m[5] = hl32_to_64(l32_from_64(pre[5]), h32_from_64(pre[4]));
    m[6] = hl32_to_64(l32_from_64(pre[6]), h32_from_64(pre[5]));
    m[7] = hl32_to_64(l32_from_64(pre[7]), h32_from_64(pre[6]));
    m[8] = hl32_to_64(                  0, h32_from_64(pre[7]));
    m[9] = hl32_to_64(                  0, swap32(lane));

    blake2b_transform(h, t, f, m, v, ARGON2_PREHASH_SEED_LENGTH + 4, BLAKE2B_FINAL);
    
    tmps[gid].blocks[lane][0] = h[0];
    tmps[gid].blocks[lane][1] = h[1];
    tmps[gid].blocks[lane][2] = h[2];
    tmps[gid].blocks[lane][3] = h[3];

    m[8] = 0;
    m[9] = 0;

    for (u32 i = 1; i < 30; i++)
    {
      m[0] = h[0];
      m[1] = h[1];
      m[2] = h[2];
      m[3] = h[3];
      m[4] = h[4];
      m[5] = h[5];
      m[6] = h[6];
      m[7] = h[7];

      h[0]  = BLAKE2B_IV_00;
      h[0] ^= BLAKE2B_BLOCK_SIZE;
      h[0] ^= BLAKE2B_FANOUT;
      h[0] ^= BLAKE2B_DEPTH;
      h[1]  = BLAKE2B_IV_01;
      h[2]  = BLAKE2B_IV_02;
      h[3]  = BLAKE2B_IV_03;
      h[4]  = BLAKE2B_IV_04;
      h[5]  = BLAKE2B_IV_05;
      h[6]  = BLAKE2B_IV_06;
      h[7]  = BLAKE2B_IV_07;

      t[0] = 0;
      t[1] = 0;
      f[0] = 0;
      f[1] = 0;

      blake2b_transform(h, t, f, m, v, BLAKE2B_BLOCK_SIZE, BLAKE2B_FINAL);

      tmps[gid].blocks[lane][i * 4 + 0] = h[0];
      tmps[gid].blocks[lane][i * 4 + 1] = h[1];
      tmps[gid].blocks[lane][i * 4 + 2] = h[2];
      tmps[gid].blocks[lane][i * 4 + 3] = h[3];
    }

    m[0] = h[0];
    m[1] = h[1];
    m[2] = h[2];
    m[3] = h[3];
    m[4] = h[4];
    m[5] = h[5];
    m[6] = h[6];
    m[7] = h[7];

    h[0]  = BLAKE2B_IV_00;
    h[0] ^= BLAKE2B_BLOCK_SIZE;
    h[0] ^= BLAKE2B_FANOUT;
    h[0] ^= BLAKE2B_DEPTH;
    h[1]  = BLAKE2B_IV_01;
    h[2]  = BLAKE2B_IV_02;
    h[3]  = BLAKE2B_IV_03;
    h[4]  = BLAKE2B_IV_04;
    h[5]  = BLAKE2B_IV_05;
    h[6]  = BLAKE2B_IV_06;
    h[7]  = BLAKE2B_IV_07;

    t[0] = 0;
    t[1] = 0;
    f[0] = 0;
    f[1] = 0;

    blake2b_transform(h, t, f, m, v, BLAKE2B_BLOCK_SIZE, BLAKE2B_FINAL);

    tmps[gid].blocks[lane][120] = h[0];
    tmps[gid].blocks[lane][121] = h[1];
    tmps[gid].blocks[lane][122] = h[2];
    tmps[gid].blocks[lane][123] = h[3];
    tmps[gid].blocks[lane][124] = h[4];
    tmps[gid].blocks[lane][125] = h[5];
    tmps[gid].blocks[lane][126] = h[6];
    tmps[gid].blocks[lane][127] = h[7];

    h[0]  = BLAKE2B_IV_00;
    h[0] ^= BLAKE2B_BLOCK_SIZE;
    h[0] ^= BLAKE2B_FANOUT;
    h[0] ^= BLAKE2B_DEPTH;
    h[1]  = BLAKE2B_IV_01;
    h[2]  = BLAKE2B_IV_02;
    h[3]  = BLAKE2B_IV_03;
    h[4]  = BLAKE2B_IV_04;
    h[5]  = BLAKE2B_IV_05;
    h[6]  = BLAKE2B_IV_06;
    h[7]  = BLAKE2B_IV_07;

    f[0] = 0;
    f[1] = 0;
    t[0] = 0;
    t[1] = 0;

    m[0] = hl32_to_64(l32_from_64(pre[0]), ARGON2_BLOCK_SIZE);
    m[1] = hl32_to_64(l32_from_64(pre[1]), h32_from_64(pre[0]));
    m[2] = hl32_to_64(l32_from_64(pre[2]), h32_from_64(pre[1]));
    m[3] = hl32_to_64(l32_from_64(pre[3]), h32_from_64(pre[2]));
    m[4] = hl32_to_64(l32_from_64(pre[4]), h32_from_64(pre[3]));
    m[5] = hl32_to_64(l32_from_64(pre[5]), h32_from_64(pre[4]));
    m[6] = hl32_to_64(l32_from_64(pre[6]), h32_from_64(pre[5]));
    m[7] = hl32_to_64(l32_from_64(pre[7]), h32_from_64(pre[6]));
    m[8] = hl32_to_64(                  1, h32_from_64(pre[7]));
    m[9] = hl32_to_64(                  0, swap32(lane));

    blake2b_transform(h, t, f, m, v, ARGON2_PREHASH_SEED_LENGTH + 4, BLAKE2B_FINAL);

    tmps[gid].blocks[lane][128] = h[0];
    tmps[gid].blocks[lane][129] = h[1];
    tmps[gid].blocks[lane][130] = h[2];
    tmps[gid].blocks[lane][131] = h[3];

    m[8] = 0;
    m[9] = 0;

    for (u32 i = 1; i < 30; i++)
    {
      m[0] = h[0];
      m[1] = h[1];
      m[2] = h[2];
      m[3] = h[3];
      m[4] = h[4];
      m[5] = h[5];
      m[6] = h[6];
      m[7] = h[7];

      h[0]  = BLAKE2B_IV_00;
      h[0] ^= BLAKE2B_BLOCK_SIZE;
      h[0] ^= BLAKE2B_FANOUT;
      h[0] ^= BLAKE2B_DEPTH;
      h[1]  = BLAKE2B_IV_01;
      h[2]  = BLAKE2B_IV_02;
      h[3]  = BLAKE2B_IV_03;
      h[4]  = BLAKE2B_IV_04;
      h[5]  = BLAKE2B_IV_05;
      h[6]  = BLAKE2B_IV_06;
      h[7]  = BLAKE2B_IV_07;

      t[0] = 0;
      t[1] = 0;
      f[0] = 0;
      f[1] = 0;

      blake2b_transform(h, t, f, m, v, BLAKE2B_BLOCK_SIZE, BLAKE2B_FINAL);

      tmps[gid].blocks[lane][i * 4 + 128] = h[0];
      tmps[gid].blocks[lane][i * 4 + 129] = h[1];
      tmps[gid].blocks[lane][i * 4 + 130] = h[2];
      tmps[gid].blocks[lane][i * 4 + 131] = h[3];
    }

    m[0] = h[0];
    m[1] = h[1];
    m[2] = h[2];
    m[3] = h[3];
    m[4] = h[4];
    m[5] = h[5];
    m[6] = h[6];
    m[7] = h[7];

    h[0]  = BLAKE2B_IV_00;
    h[0] ^= BLAKE2B_BLOCK_SIZE;
    h[0] ^= BLAKE2B_FANOUT;
    h[0] ^= BLAKE2B_DEPTH;
    h[1]  = BLAKE2B_IV_01;
    h[2]  = BLAKE2B_IV_02;
    h[3]  = BLAKE2B_IV_03;
    h[4]  = BLAKE2B_IV_04;
    h[5]  = BLAKE2B_IV_05;
    h[6]  = BLAKE2B_IV_06;
    h[7]  = BLAKE2B_IV_07;

    t[0] = 0;
    t[1] = 0;
    f[0] = 0;
    f[1] = 0;

    blake2b_transform(h, t, f, m, v, BLAKE2B_BLOCK_SIZE, BLAKE2B_FINAL);

    tmps[gid].blocks[lane][248] = h[0];
    tmps[gid].blocks[lane][249] = h[1];
    tmps[gid].blocks[lane][250] = h[2];
    tmps[gid].blocks[lane][251] = h[3];
    tmps[gid].blocks[lane][252] = h[4];
    tmps[gid].blocks[lane][253] = h[5];
    tmps[gid].blocks[lane][254] = h[6];
    tmps[gid].blocks[lane][255] = h[7];
  }

  
  u64 tmp[1][256];

  for (u32 i = 0; i < 256; i++) tmp[0][i] = tmps[gid].blocks[0][i];

  for (u32 i = 0; i < 1; i++)
  {
    printf("\n");
    for (u32 j = 0; j < 1024; j++)
    { 
      printf("%02x", ((const u8 *)tmp[i])[j]); 
    }

  }

  printf("\n\n");
  printf("last 32... %08x%08x%08x%08x\n", ((const u32 *)&tmp[0][248])[0], ((const u32 *)&tmp[0][248])[1], ((const u32 *)&tmp[0][248])[2], ((const u32 *)&tmp[0][248])[3] );
}

__kernel void m15500_loop (__global pw_t *pws, __global const kernel_rule_t *rules_buf, __global const comb_t *combs_buf, __global const bf_t *bfs_buf, __global argon2_tmp_t *tmps, __global void *hooks, __global const u32 *bitmaps_buf_s1_a, __global const u32 *bitmaps_buf_s1_b, __global const u32 *bitmaps_buf_s1_c, __global const u32 *bitmaps_buf_s1_d, __global const u32 *bitmaps_buf_s2_a, __global const u32 *bitmaps_buf_s2_b, __global const u32 *bitmaps_buf_s2_c, __global const u32 *bitmaps_buf_s2_d, __global plain_t *plains_buf, __global const digest_t *digests_buf, __global u32 *hashes_shown, __global const salt_t *salt_bufs, __global argon2_t *esalt_bufs, __global u32 *d_return_buf, __global u32 *d_scryptV0_buf, __global u32 *d_scryptV1_buf, __global u32 *d_scryptV2_buf, __global u32 *d_scryptV3_buf, const u32 bitmap_mask, const u32 bitmap_shift1, const u32 bitmap_shift2, const u32 salt_pos, const u32 loop_pos, const u32 loop_cnt, const u32 rules_cnt, const u32 digests_cnt, const u32 digests_offset, const u32 combs_mode, const u32 gid_max)
{
  /**
   * base
   */

  const u32 gid = get_global_id (0);

  if ((gid * VECT_SIZE) >= gid_max) return;


  /**
   * iter1
   */



  /**
   * Update tmps
   */

}

__kernel void m15500_comp (__global pw_t *pws, __global const kernel_rule_t *rules_buf, __global const comb_t *combs_buf, __global const bf_t *bfs_buf, __global argon2_tmp_t *tmps, __global void *hooks, __global const u32 *bitmaps_buf_s1_a, __global const u32 *bitmaps_buf_s1_b, __global const u32 *bitmaps_buf_s1_c, __global const u32 *bitmaps_buf_s1_d, __global const u32 *bitmaps_buf_s2_a, __global const u32 *bitmaps_buf_s2_b, __global const u32 *bitmaps_buf_s2_c, __global const u32 *bitmaps_buf_s2_d, __global plain_t *plains_buf, __global const digest_t *digests_buf, __global u32 *hashes_shown, __global const salt_t *salt_bufs, __global argon2_t *esalt_bufs, __global u32 *d_return_buf, __global u32 *d_scryptV0_buf, __global u32 *d_scryptV1_buf, __global u32 *d_scryptV2_buf, __global u32 *d_scryptV3_buf, const u32 bitmap_mask, const u32 bitmap_shift1, const u32 bitmap_shift2, const u32 salt_pos, const u32 loop_pos, const u32 loop_cnt, const u32 rules_cnt, const u32 digests_cnt, const u32 digests_offset, const u32 combs_mode, const u32 gid_max)
{
  const u32 gid = get_global_id (0);
  const u32 lid = get_local_id (0);
  const u32 lsz = get_local_size (0);

  /**
   * Mark hashes
   */

}

