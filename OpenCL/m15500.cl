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

#define BLAKE2B_FINAL   1
#define BLAKE2B_UPDATE  0

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

  for (u8 i = 0; i < 8; i++)
  {
    printf("h1[%d] = %08x%08x\n", i, h32_from_64(h[i]), l32_from_64(h[i]));
  }

  printf("t %d, f, %d\n", t[0], f[0]);

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

  for (u8 i = 0; i < 8; i++)
  {
    printf("h2[%d] = %08x%08x\n", i, h32_from_64(h[i]), l32_from_64(h[i]));
  }

}



__kernel void m15500_init (__global pw_t *pws, __global const kernel_rule_t *rules_buf, __global const comb_t *combs_buf, __global const bf_t *bfs_buf, __global argon2_tmp_t *tmps, __global void *hooks, __global const u32 *bitmaps_buf_s1_a, __global const u32 *bitmaps_buf_s1_b, __global const u32 *bitmaps_buf_s1_c, __global const u32 *bitmaps_buf_s1_d, __global const u32 *bitmaps_buf_s2_a, __global const u32 *bitmaps_buf_s2_b, __global const u32 *bitmaps_buf_s2_c, __global const u32 *bitmaps_buf_s2_d, __global plain_t *plains_buf, __global const digest_t *digests_buf, __global u32 *hashes_shown, __global const salt_t *salt_bufs, __global argon2_t *esalt_bufs, __global u32 *d_return_buf, __global u32 *d_scryptV0_buf, __global u32 *d_scryptV1_buf, __global u32 *d_scryptV2_buf, __global u32 *d_scryptV3_buf, const u32 bitmap_mask, const u32 bitmap_shift1, const u32 bitmap_shift2, const u32 salt_pos, const u32 loop_pos, const u32 loop_cnt, const u32 rules_cnt, const u32 digests_cnt, const u32 digests_offset, const u32 combs_mode, const u32 gid_max)
{
  /**
   * base
   */

  const u32 gid   = get_global_id (0);

  salt_t   salt   =  salt_bufs[salt_pos];
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
   * Init Blake2
   */

  u64x m[16] = { 0 };
  u64x v[16] = { 0 };
  u64x h[ 8] = { 0 };
  u64x f[ 2] = { 0 };
  u64x t[ 2] = { 0 };

/*
  m[ 0] = swap32( argon2.p );
  m[ 1] = swap32( argon2.l );   
  m[ 2] = swap32( argon2.m );
  m[ 3] = swap32( argon2.t );
  m[ 4] = swap32( argon2.v );
  m[ 5] = swap32( argon2.y );
  m[ 6] = swap32( pw_len   );
  m[ 7] = swap32( w0[0]    );
  m[ 8] = swap32( w0[1]    );
  m[ 9] = swap32( w0[2]    );
  m[10] = swap32( w0[3]    );
  m[11] = swap32( w0[4]    );
  m[12] = swap32( w0[5]    );
  m[13] = swap32( w0[6]    );
  m[14] = swap32( w0[7]    );
*/

  u8x p[128] = { 0 };
  u8  index  = 0; 

  p[ index++ ] = argon2.p >>  0 & 0xff;
  p[ index++ ] = argon2.p >>  8 & 0xff;
  p[ index++ ] = argon2.p >> 16 & 0xff;
  p[ index++ ] = argon2.p >> 24 & 0xff;
  p[ index++ ] = argon2.l >>  0 & 0xff;
  p[ index++ ] = argon2.l >>  8 & 0xff;
  p[ index++ ] = argon2.l >> 16 & 0xff;
  p[ index++ ] = argon2.l >> 24 & 0xff;
  p[ index++ ] = argon2.m >>  0 & 0xff;
  p[ index++ ] = argon2.m >>  8 & 0xff;
  p[ index++ ] = argon2.m >> 16 & 0xff;
  p[ index++ ] = argon2.m >> 24 & 0xff;
  p[ index++ ] = argon2.t >>  0 & 0xff;
  p[ index++ ] = argon2.t >>  8 & 0xff;
  p[ index++ ] = argon2.t >> 16 & 0xff;
  p[ index++ ] = argon2.t >> 24 & 0xff;
  p[ index++ ] = argon2.v >>  0 & 0xff;
  p[ index++ ] = argon2.v >>  8 & 0xff;
  p[ index++ ] = argon2.v >> 16 & 0xff;
  p[ index++ ] = argon2.v >> 24 & 0xff;
  p[ index++ ] = argon2.y >>  0 & 0xff;
  p[ index++ ] = argon2.y >>  8 & 0xff;
  p[ index++ ] = argon2.y >> 16 & 0xff;
  p[ index++ ] = argon2.y >> 24 & 0xff;
  p[ index++ ] = pw_len   >>  0 & 0xff;
  p[ index++ ] = pw_len   >>  8 & 0xff;
  p[ index++ ] = pw_len   >> 16 & 0xff;
  p[ index++ ] = pw_len   >> 24 & 0xff;
  
  for (u8 i = 0; i < pw_len; i++)
  {
    p[ index++ ] = w0[i / 4] >> ((i % 4) * 8) & 0xff;
  }

  p[ index++ ] = salt.salt_len >>  0 & 0xff;
  p[ index++ ] = salt.salt_len >>  8 & 0xff;
  p[ index++ ] = salt.salt_len >> 16 & 0xff;
  p[ index++ ] = salt.salt_len >> 24 & 0xff;

  for (u8 i = 0; i < salt.salt_len; i++)
  {  
    p[ index++ ] = salt.salt_buf[i / 4] >> ((i % 4) * 8) & 0xff;
  }

  for (u8 i = 0; i < 16; i++)
  {
    m[i] =  p[i * 8 + 4] <<  0
         |  p[i * 8 + 5] <<  8
         |  p[i * 8 + 6] << 16
         |  p[i * 8 + 7] << 24;

    m[i] =  m[i] << 32;

    m[i] |= p[i * 8 + 0] <<  0
         |  p[i * 8 + 1] <<  8
         |  p[i * 8 + 2] << 16
         |  p[i * 8 + 3] << 24;
  }
 
  h[0]  = BLAKE2B_IV_00; 
  h[0] ^= 0x00000040;      // outlen 
  h[0] ^= 0x00010000;      // fanout
  h[0] ^= 0x01000000;      // depth
  h[1]  = BLAKE2B_IV_01;
  h[2]  = BLAKE2B_IV_02;
  h[3]  = BLAKE2B_IV_03;
  h[4]  = BLAKE2B_IV_04;
  h[5]  = BLAKE2B_IV_05;
  h[6]  = BLAKE2B_IV_06;
  h[7]  = BLAKE2B_IV_07;

  u64x total_size = (10 * 4) + pw_len + salt.salt_len;

  blake2b_transform(h, t, f, m, v, total_size, BLAKE2B_FINAL);
  

  /**
   * Store tmps
   */
  
  printf("\n");
  printf("Digest:       %08x%08x%08x%08x\n", digests_buf[digests_offset].digest_buf[0], digests_buf[digests_offset].digest_buf[1], digests_buf[digests_offset].digest_buf[2], digests_buf[digests_offset].digest_buf[3]);
  printf("Password:     %08x%08x%08x%08x\n", w0[0], w0[1], w0[2], w0[3]);
  printf("Passw Length: %d\n", pw_len);
  printf("Salt:         %08x%08x, salt_len: %d\n", salt.salt_buf[0], salt_bufs[salt_pos].salt_buf[1], salt_bufs[salt_pos].salt_len);
  printf("Esalt:        y=%d, v=%d, m=%d, t=%d, p=%d, l=%d\n", esalt_bufs[salt_pos].y, esalt_bufs[salt_pos].v, esalt_bufs[salt_pos].m, esalt_bufs[salt_pos].t, esalt_bufs[salt_pos].p, esalt_bufs[salt_pos].l);

  u8 *m_ptr = (const u8 *)m;

  //printf("\nm[x]:\n");
  for (u8 i = 0; i < 64; i++) printf("%02x", m_ptr[i]);
  printf("\n01000000200000000010000003000000130000000000000005000000706173737708000000313233343536373800000000000000000000000000000000000000\n");  
  printf("\nHash: ");
  for (u8 i = 0; i < 16; i++) printf("%08x", ((const u32 *) h)[i]);
  printf("\n");
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

