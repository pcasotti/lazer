/*
 * \file   abdlop_arith.c
 * \brief  abdlop pipeline over the pluggable arith interface.
 *
 * Self-contained: depends only on arith.h and standard/crypto libs (no
 * lazer.h).  Mirrors src/abdlop.c; ring multiplications go through arith
 * ops, while sampling, encoding, hashing and rejection are implemented
 * locally with deterministic, bit-compatible algorithms.
 */
#include "abdlop_arith.h"

#include <endian.h>
#include <immintrin.h>
#include <math.h>
#include <mpfr.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ========================= coeff helpers ========================= */

static long long
_canonq (long long v, unsigned long long q)
{
  v %= (long long)q;
  if (v < 0) v += (long long)q;
  return v;
}

static long long
_balq (long long c, unsigned long long q)
{
  /* c expected canonical [0,q) */
  long long v = c % (long long)q;
  if (v < 0) v += (long long)q;
  const long long pivot = (long long)((q - 1) >> 1);
  return v > pivot ? v - (long long)q : v;
}

static void
_decomp_c (long long c, long long *r1p, long long *r0p,
           const abdlop_arith_params_t *ap)
{
  const long long gamma = ap->gamma, gb2 = ap->gammaby2;
  long long r1;
  c = _canonq (c, ap->q);
  long long r0 = c % gamma;
  if (r0 > gb2) r0 -= gamma;
  r1 = c - r0;
  if (r1 == ap->qminus1) { r1 = 0; r0 -= 1; }
  else r1 /= gamma;
  *r1p = r1; *r0p = r0;
}

static long long
_p2r_c (long long c, long long *t2p, const abdlop_arith_params_t *ap)
{
  const long long pd = ap->pow2D, pd2 = ap->pow2Dby2;
  long long r0, t1, t2;
  c = _canonq (c, ap->q);
  r0 = c & (pd - 1);
  if (r0 > pd2) r0 -= pd;
  t1 = (c - r0) >> ap->D;
  t2 = c - (t1 << ap->D);
  *t2p = _canonq (t2, ap->q);
  return t1;
}

static long long
_use_ghint_c (long long r, long long y, const abdlop_arith_params_t *ap)
{
  long long r1, r0, s;
  _decomp_c (r, &r1, &r0, ap); (void)r0;
  s = r1 + y;
  s %= ap->m;
  if (s < 0) s += ap->m;
  return s;
}

static long long
_make_ghint_c (long long r, long long z, const abdlop_arith_params_t *ap)
{
  const long long m = ap->m, mb2 = ap->mby2;
  long long s, r1a, r0a, r1b, r0b, h;
  _decomp_c (r, &r1a, &r0a, ap); (void)r0a;
  s = _canonq (r + z, ap->q);
  _decomp_c (s, &r1b, &r0b, ap); (void)r0b;
  h = r1b - r1a;
  h %= m; if (h < 0) h += m;
  if (h < 0) { if (-h >= mb2) h += m; }
  else { if (h > mb2) h -= m; }
  return h;
}

static __int128
_l2sqr (abdlop_arith_engine_t *e, const arith_vec_ptr v, unsigned int n,
        const abdlop_arith_params_t *ap)
{
  __int128 s = 0; unsigned int i, k;
  for (i = 0; i < n; i++)
    for (k = 0; k < ap->d; k++) {
      long long b = _balq (e->ops->vec_get_coeff (v, i, k), ap->q);
      s += (__int128)b * b;
    }
  return s;
}

static long long
_linf (abdlop_arith_engine_t *e, const arith_vec_ptr v, unsigned int n,
       const abdlop_arith_params_t *ap)
{
  long long mx = 0; unsigned int i, k;
  for (i = 0; i < n; i++)
    for (k = 0; k < ap->d; k++) {
      long long b = _balq (e->ops->vec_get_coeff (v, i, k), ap->q);
      if (b < 0) b = -b;
      if (b > mx) mx = b;
    }
  return mx;
}

static void
_v2bal (abdlop_arith_engine_t *e, const abdlop_arith_params_t *ap,
        long long *out, const arith_vec_ptr v, unsigned int n)
{
  unsigned int i, k;
  for (i = 0; i < n; i++)
    for (k = 0; k < ap->d; k++)
      out[i * ap->d + k] = _balq (e->ops->vec_get_coeff (v, i, k), ap->q);
}

/* ========================= local SHAKE128 ========================= */
#define AA_SHAKE128_RATE 168
#define AA_NROUNDS 24
#define AA_ROL(a, offset) (((a) << (offset)) ^ ((a) >> (64 - (offset))))

typedef struct { uint64_t s[25]; unsigned int pos; int final; } aa_shake128_state;

static uint64_t aa_load64 (const uint8_t x[8]) { unsigned int i; uint64_t r=0; for(i=0;i<8;i++) r|=(uint64_t)x[i]<<8*i; return r; }
static void aa_store64 (uint8_t x[8], uint64_t u) { unsigned int i; for(i=0;i<8;i++) x[i]=(uint8_t)(u>>8*i); }
static const uint64_t aa_KeccakF_RoundConstants[AA_NROUNDS] = {
  0x0000000000000001ULL,0x0000000000008082ULL,0x800000000000808aULL,0x8000000080008000ULL,
  0x000000000000808bULL,0x0000000080000001ULL,0x8000000080008081ULL,0x8000000000008009ULL,
  0x000000000000008aULL,0x0000000000000088ULL,0x0000000080008009ULL,0x000000008000000aULL,
  0x000000008000808bULL,0x800000000000008bULL,0x8000000000008089ULL,0x8000000000008003ULL,
  0x8000000000008002ULL,0x8000000000000080ULL,0x000000000000800aULL,0x800000008000000aULL,
  0x8000000080008081ULL,0x8000000000008080ULL,0x0000000080000001ULL,0x8000000080008008ULL };

static void
aa_KeccakF1600_StatePermute (uint64_t state[25])
{
  int round;
  uint64_t Aba,Abe,Abi,Abo,Abu,Aga,Age,Agi,Ago,Agu,Aka,Ake,Aki,Ako,Aku,Ama,Ame,Ami,Amo,Amu,Asa,Ase,Asi,Aso,Asu;
  uint64_t BCa,BCe,BCi,BCo,BCu,Da,De,Di,Do,Du;
  uint64_t Eba,Ebe,Ebi,Ebo,Ebu,Ega,Ege,Egi,Ego,Egu,Eka,Eke,Eki,Eko,Eku,Ema,Eme,Emi,Emo,Emu,Esa,Ese,Esi,Eso,Esu;
  Aba=state[0];Abe=state[1];Abi=state[2];Abo=state[3];Abu=state[4];Aga=state[5];Age=state[6];Agi=state[7];Ago=state[8];Agu=state[9];
  Aka=state[10];Ake=state[11];Aki=state[12];Ako=state[13];Aku=state[14];Ama=state[15];Ame=state[16];Ami=state[17];Amo=state[18];Amu=state[19];
  Asa=state[20];Ase=state[21];Asi=state[22];Aso=state[23];Asu=state[24];
  for(round=0;round<AA_NROUNDS;round+=2){
    BCa=Aba^Aga^Aka^Ama^Asa; BCe=Abe^Age^Ake^Ame^Ase; BCi=Abi^Agi^Aki^Ami^Asi; BCo=Abo^Ago^Ako^Amo^Aso; BCu=Abu^Agu^Aku^Amu^Asu;
    Da=BCu^AA_ROL(BCe,1); De=BCa^AA_ROL(BCi,1); Di=BCe^AA_ROL(BCo,1); Do=BCi^AA_ROL(BCu,1); Du=BCo^AA_ROL(BCa,1);
    Aba^=Da;BCa=Aba;Age^=De;BCe=AA_ROL(Age,44);Aki^=Di;BCi=AA_ROL(Aki,43);Amo^=Do;BCo=AA_ROL(Amo,21);Asu^=Du;BCu=AA_ROL(Asu,14);
    Eba=BCa^((~BCe)&BCi);Eba^=aa_KeccakF_RoundConstants[round];Ebe=BCe^((~BCi)&BCo);Ebi=BCi^((~BCo)&BCu);Ebo=BCo^((~BCu)&BCa);Ebu=BCu^((~BCa)&BCe);
    Abo^=Do;BCa=AA_ROL(Abo,28);Agu^=Du;BCe=AA_ROL(Agu,20);Aka^=Da;BCi=AA_ROL(Aka,3);Ame^=De;BCo=AA_ROL(Ame,45);Asi^=Di;BCu=AA_ROL(Asi,61);
    Ega=BCa^((~BCe)&BCi);Ege=BCe^((~BCi)&BCo);Egi=BCi^((~BCo)&BCu);Ego=BCo^((~BCu)&BCa);Egu=BCu^((~BCa)&BCe);
    Abe^=De;BCa=AA_ROL(Abe,1);Agi^=Di;BCe=AA_ROL(Agi,6);Ako^=Do;BCi=AA_ROL(Ako,25);Amu^=Du;BCo=AA_ROL(Amu,8);Asa^=Da;BCu=AA_ROL(Asa,18);
    Eka=BCa^((~BCe)&BCi);Eke=BCe^((~BCi)&BCo);Eki=BCi^((~BCo)&BCu);Eko=BCo^((~BCu)&BCa);Eku=BCu^((~BCa)&BCe);
    Abu^=Du;BCa=AA_ROL(Abu,27);Aga^=Da;BCe=AA_ROL(Aga,36);Ake^=De;BCi=AA_ROL(Ake,10);Ami^=Di;BCo=AA_ROL(Ami,15);Aso^=Do;BCu=AA_ROL(Aso,56);
    Ema=BCa^((~BCe)&BCi);Eme=BCe^((~BCi)&BCo);Emi=BCi^((~BCo)&BCu);Emo=BCo^((~BCu)&BCa);Emu=BCu^((~BCa)&BCe);
    Abi^=Di;BCa=AA_ROL(Abi,62);Ago^=Do;BCe=AA_ROL(Ago,55);Aku^=Du;BCi=AA_ROL(Aku,39);Ama^=Da;BCo=AA_ROL(Ama,41);Ase^=De;BCu=AA_ROL(Ase,2);
    Esa=BCa^((~BCe)&BCi);Ese=BCe^((~BCi)&BCo);Esi=BCi^((~BCo)&BCu);Eso=BCo^((~BCu)&BCa);Esu=BCu^((~BCa)&BCe);
    BCa=Eba^Ega^Eka^Ema^Esa; BCe=Ebe^Ege^Eke^Eme^Ese; BCi=Ebi^Egi^Eki^Emi^Esi; BCo=Ebo^Ego^Eko^Emo^Eso; BCu=Ebu^Egu^Eku^Emu^Esu;
    Da=BCu^AA_ROL(BCe,1); De=BCa^AA_ROL(BCi,1); Di=BCe^AA_ROL(BCo,1); Do=BCi^AA_ROL(BCu,1); Du=BCo^AA_ROL(BCa,1);
    Eba^=Da;BCa=Eba;Ege^=De;BCe=AA_ROL(Ege,44);Eki^=Di;BCi=AA_ROL(Eki,43);Emo^=Do;BCo=AA_ROL(Emo,21);Esu^=Du;BCu=AA_ROL(Esu,14);
    Aba=BCa^((~BCe)&BCi);Aba^=aa_KeccakF_RoundConstants[round+1];Abe=BCe^((~BCi)&BCo);Abi=BCi^((~BCo)&BCu);Abo=BCo^((~BCu)&BCa);Abu=BCu^((~BCa)&BCe);
    Ebo^=Do;BCa=AA_ROL(Ebo,28);Egu^=Du;BCe=AA_ROL(Egu,20);Eka^=Da;BCi=AA_ROL(Eka,3);Eme^=De;BCo=AA_ROL(Eme,45);Esi^=Di;BCu=AA_ROL(Esi,61);
    Aga=BCa^((~BCe)&BCi);Age=BCe^((~BCi)&BCo);Agi=BCi^((~BCo)&BCu);Ago=BCo^((~BCu)&BCa);Agu=BCu^((~BCa)&BCe);
    Ebe^=De;BCa=AA_ROL(Ebe,1);Egi^=Di;BCe=AA_ROL(Egi,6);Eko^=Do;BCi=AA_ROL(Eko,25);Emu^=Du;BCo=AA_ROL(Emu,8);Esa^=Da;BCu=AA_ROL(Esa,18);
    Aka=BCa^((~BCe)&BCi);Ake=BCe^((~BCi)&BCo);Aki=BCi^((~BCo)&BCu);Ako=BCo^((~BCu)&BCa);Aku=BCu^((~BCa)&BCe);
    Ebu^=Du;BCa=AA_ROL(Ebu,27);Ega^=Da;BCe=AA_ROL(Ega,36);Eke^=De;BCi=AA_ROL(Eke,10);Emi^=Di;BCo=AA_ROL(Emi,15);Eso^=Do;BCu=AA_ROL(Eso,56);
    Ama=BCa^((~BCe)&BCi);Ame=BCe^((~BCi)&BCo);Ami=BCi^((~BCo)&BCu);Amo=BCo^((~BCu)&BCa);Amu=BCu^((~BCa)&BCe);
    Ebi^=Di;BCa=AA_ROL(Ebi,62);Ego^=Do;BCe=AA_ROL(Ego,55);Eku^=Du;BCi=AA_ROL(Eku,39);Ema^=Da;BCo=AA_ROL(Ema,41);Ese^=De;BCu=AA_ROL(Ese,2);
    Asa=BCa^((~BCe)&BCi);Ase=BCe^((~BCi)&BCo);Asi=BCi^((~BCo)&BCu);Aso=BCo^((~BCu)&BCa);Asu=BCu^((~BCa)&BCe);
  }
  state[0]=Aba;state[1]=Abe;state[2]=Abi;state[3]=Abo;state[4]=Abu;state[5]=Aga;state[6]=Age;state[7]=Agi;state[8]=Ago;state[9]=Agu;
  state[10]=Aka;state[11]=Ake;state[12]=Aki;state[13]=Ako;state[14]=Aku;state[15]=Ama;state[16]=Ame;state[17]=Ami;state[18]=Amo;state[19]=Amu;
  state[20]=Asa;state[21]=Ase;state[22]=Asi;state[23]=Aso;state[24]=Asu;
}
static void aa_keccak_init(aa_shake128_state *st){ unsigned int i; for(i=0;i<25;i++) st->s[i]=0; st->pos=0; st->final=0; }
static unsigned int aa_keccak_absorb(uint64_t s[25],unsigned int r,unsigned int pos,const uint8_t *m,size_t mlen){
  unsigned int i; uint8_t t[8]={0};
  if(pos&7){ i=pos&7; while(i<8&&mlen>0){ t[i++]=*m++; mlen--; pos++; } s[(pos-i)/8]^=aa_load64(t); }
  if(pos&&mlen>=r-pos){ for(i=0;i<(r-pos)/8;i++) s[pos/8+i]^=aa_load64(m+8*i); m+=r-pos; mlen-=r-pos; pos=0; aa_KeccakF1600_StatePermute(s); }
  while(mlen>=r){ for(i=0;i<r/8;i++) s[i]^=aa_load64(m+8*i); m+=r; mlen-=r; aa_KeccakF1600_StatePermute(s); }
  for(i=0;i<mlen/8;i++) s[pos/8+i]^=aa_load64(m+8*i); m+=8*i; mlen-=8*i; pos+=8*i;
  if(mlen){ for(i=0;i<8;i++) t[i]=0; for(i=0;i<(unsigned int)mlen;i++) t[i]=m[i]; s[pos/8]^=aa_load64(t); pos+=(unsigned int)mlen; }
  return pos;
}
static void aa_keccak_finalize(uint64_t s[25],unsigned int r,unsigned int pos,uint8_t p){ unsigned int i=pos>>3,j=pos&7; s[i]^=(uint64_t)p<<8*j; s[r/8-1]^=1ULL<<63; }
static unsigned int aa_keccak_squeeze(uint8_t *out,size_t outlen,uint64_t s[25],unsigned int r,unsigned int pos){
  unsigned int i; uint8_t t[8];
  if(pos&7){ aa_store64(t,s[pos/8]); i=pos&7; while(i<8&&outlen>0){ *out++=t[i++]; outlen--; pos++; } }
  if(pos&&outlen>=r-pos){ for(i=0;i<(r-pos)/8;i++) aa_store64(out+8*i,s[pos/8+i]); out+=r-pos; outlen-=r-pos; pos=0; }
  while(outlen>=r){ aa_KeccakF1600_StatePermute(s); for(i=0;i<r/8;i++) aa_store64(out+8*i,s[i]); out+=r; outlen-=r; }
  if(!outlen) return pos; else if(!pos) aa_KeccakF1600_StatePermute(s);
  for(i=0;i<outlen/8;i++) aa_store64(out+8*i,s[pos/8+i]); out+=8*i; outlen-=8*i; pos+=8*i;
  aa_store64(t,s[pos/8]); for(i=0;i<outlen;i++) out[i]=t[i]; pos+=(unsigned int)outlen; return pos;
}
static void aa_shake128_init(aa_shake128_state *st){ aa_keccak_init(st); }
static void aa_shake128_absorb(aa_shake128_state *st,const uint8_t *in,size_t len){ st->pos=aa_keccak_absorb(st->s,AA_SHAKE128_RATE,st->pos,in,len); }
static void aa_shake128_squeeze(aa_shake128_state *st,uint8_t *out,size_t len){
  if(st->final==0){ aa_keccak_finalize(st->s,AA_SHAKE128_RATE,st->pos,0x1F); st->pos=0; st->final=1; }
  st->pos=aa_keccak_squeeze(out,len,st->s,AA_SHAKE128_RATE,st->pos);
}
static void aa_shake128_clear(aa_shake128_state *st){ memset(st,0,sizeof(*st)); }

/* ========================= local AES256CTR (amd64/AES-NI) + RNG ========================= */
typedef struct {
  __m128i rkeys[16];
  __attribute__((aligned(16))) uint8_t n2[16];
  __attribute__((aligned(16))) uint8_t cache[8 * 16];
  uint8_t *cache_ptr;
  unsigned int nbytes;
} aa_aes256ctr_state;

static inline void aa_aesni_key256_expand(const unsigned char *key, __m128i rkeys[16]){
  __m128i key0=_mm_loadu_si128((const __m128i*)(key+0));
  __m128i key1=_mm_loadu_si128((const __m128i*)(key+16));
  __m128i temp0,temp1,temp2,temp4; int idx=0;
  rkeys[idx++]=key0; temp0=key0; temp2=key1;
#define AA_BLOCK1(IMM) temp1=_mm_aeskeygenassist_si128(temp2,IMM); rkeys[idx++]=temp2; temp4=_mm_slli_si128(temp0,4); temp0=_mm_xor_si128(temp0,temp4); temp4=_mm_slli_si128(temp0,8); temp0=_mm_xor_si128(temp0,temp4); temp1=_mm_shuffle_epi32(temp1,0xff); temp0=_mm_xor_si128(temp0,temp1)
#define AA_BLOCK2(IMM) temp1=_mm_aeskeygenassist_si128(temp0,IMM); rkeys[idx++]=temp0; temp4=_mm_slli_si128(temp2,4); temp2=_mm_xor_si128(temp2,temp4); temp4=_mm_slli_si128(temp2,8); temp2=_mm_xor_si128(temp2,temp4); temp1=_mm_shuffle_epi32(temp1,0xaa); temp2=_mm_xor_si128(temp2,temp1)
  AA_BLOCK1(0x01); AA_BLOCK2(0x01); AA_BLOCK1(0x02); AA_BLOCK2(0x02);
  AA_BLOCK1(0x04); AA_BLOCK2(0x04); AA_BLOCK1(0x08); AA_BLOCK2(0x08);
  AA_BLOCK1(0x10); AA_BLOCK2(0x10); AA_BLOCK1(0x20); AA_BLOCK2(0x20);
  AA_BLOCK1(0x40); rkeys[idx++]=temp0;
#undef AA_BLOCK1
#undef AA_BLOCK2
}
#define AA_NVx(a) __m128i nv##a;
#define AA_NVx_NOWRAP(a) nv##a=_mm_shuffle_epi8(_mm_add_epi64(nv0i,_mm_set_epi64x(a,0)),_mm_set_epi8(8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7))
#define AA_NVx_WRAP(a) {__m128i ad=_mm_add_epi64(nv0i,_mm_set_epi64x(a,a>=wrapnumber?1:0)); nv##a=_mm_shuffle_epi8(ad,_mm_set_epi8(8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7));}
#define AA_TEMPx(a) __m128i temp##a=_mm_xor_si128(nv##a,rkeys[0])
#define AA_AESENCx(a) temp##a=_mm_aesenc_si128(temp##a,rkeys[i]);
#define AA_AESENCLASTx(a) temp##a=_mm_aesenclast_si128(temp##a,rkeys[14]);
#define AA_STOREx(a) _mm_store_si128((__m128i*)(out+(a*16)),temp##a);
#define AA_M8(X) X(0);X(1);X(2);X(3);X(4);X(5);X(6);X(7)
static inline void aa_aesni_encrypt8(unsigned char *out,unsigned char *n,__m128i rkeys[16]){
  __m128i nv0i=_mm_load_si128((const __m128i*)n); long long nl=*(long long*)&n[8];
  AA_M8(AA_NVx);
  if((nl<0)&&(nl+8)>=0){ int wrapnumber=(int)(8-(nl+8)); AA_M8(AA_NVx_WRAP); _mm_storeu_si128((__m128i*)n,_mm_add_epi64(nv0i,_mm_set_epi64x(8,1))); }
  else { AA_M8(AA_NVx_NOWRAP); _mm_storeu_si128((__m128i*)n,_mm_add_epi64(nv0i,_mm_set_epi64x(8,0))); }
  int i; AA_M8(AA_TEMPx); for(i=1;i<14;i++){ AA_M8(AA_AESENCx); } AA_M8(AA_AESENCLASTx); AA_M8(AA_STOREx);
}
#undef AA_NVx
#undef AA_NVx_NOWRAP
#undef AA_NVx_WRAP
#undef AA_TEMPx
#undef AA_AESENCx
#undef AA_AESENCLASTx
#undef AA_STOREx
#undef AA_M8
static void aa_aes256ctr_init(aa_aes256ctr_state *state,const uint8_t key[32],const uint8_t nonce[16]){
  state->cache_ptr=NULL; state->nbytes=0; aa_aesni_key256_expand(key,state->rkeys);
  (*(unsigned long long *)&state->n2[8]) = __builtin_bswap64((*(const unsigned long long *)&nonce[8]));
  (*(unsigned long long *)&state->n2[0]) = __builtin_bswap64((*(const unsigned long long *)&nonce[0]));
}
static void aa_aes256ctr_stream(aa_aes256ctr_state *state,uint8_t *out,size_t outlen){
  size_t i,j,mj,len; len = outlen < state->nbytes ? outlen : state->nbytes;
  for(i=0;i<len;i++) out[i]=state->cache_ptr[i];
  state->cache_ptr += len; state->nbytes -= (unsigned int)len; out += len; outlen -= len;
  for(i=0;i<outlen;){ aa_aesni_encrypt8(state->cache,state->n2,state->rkeys); mj=128; if((i+mj)>=outlen){ mj=outlen-i; state->nbytes=(unsigned int)(128-mj); state->cache_ptr=state->cache+mj; } for(j=0;j<mj;j++) out[i+j]=state->cache[j]; i+=mj; if(i>=outlen) break; }
}
static void aa_aes256ctr_clear(aa_aes256ctr_state *state){ memset(state,0,sizeof(*state)); }

typedef struct { aa_aes256ctr_state st; } aa_rng_state;
static void aa_rng_init(aa_rng_state *r,const uint8_t seed[32],uint64_t dom){
  uint64_t nonce[2]; dom=htole64(dom); nonce[0]=dom; nonce[1]=0;
  aa_aes256ctr_init(&r->st, seed, (uint8_t*)nonce);
}
static void aa_rng_urandom(aa_rng_state *r,uint8_t *out,size_t outlen){ if(outlen>0) aa_aes256ctr_stream(&r->st,out,outlen); }
static void aa_rng_clear(aa_rng_state *r){ aa_aes256ctr_clear(&r->st); }

/* ========================= local coder (urandom bit-packing) ========================= */
typedef struct { uint8_t *out; unsigned int byte_off; unsigned int bit_off; } aa_coder;
static void aa_coder_begin(aa_coder *s,uint8_t *out){ s->out=out; s->byte_off=0; s->bit_off=0; }
static unsigned int aa_coder_offset(aa_coder *s){ return (s->byte_off<<3)+s->bit_off; }
static inline void aa_inc_idx_zero(uint8_t *buf,unsigned int *byte,unsigned int *bit){ (*bit)++; if(*bit>7){ *bit=0; (*byte)++; buf[*byte]=0; } }
static void aa_coder_end(aa_coder *s){ uint8_t *b=s->out; unsigned int bit=s->bit_off; b[0]&=(uint8_t)~((uint8_t)(~0)<<bit); b[0]|=(uint8_t)(1<<bit); s->bit_off=0; s->byte_off++; }
/* encode a vector: nelems polys * d coeffs, mbits each, element-then-coeff order (matches coder_enc_urandom3) */
static void
aa_encode_urandom_vec(aa_coder *s, abdlop_arith_engine_t *e,
                      const arith_vec_ptr v, unsigned int nelems,
                      unsigned int d, unsigned int mbits)
{
  unsigned int i,k,j; uint8_t *b=s->out; unsigned int bit=s->bit_off; unsigned int nbytes=0;
  b[0] &= (uint8_t)~((uint8_t)(~0) << bit);
  for(i=0;i<nelems;i++)
    for(k=0;k<d;k++){
      uint64_t val=(uint64_t)e->ops->vec_get_coeff(v,i,k);
      for(j=0;j<mbits;j++){ b[nbytes] |= (uint8_t)((val&1ULL)<<bit); val>>=1; aa_inc_idx_zero(b,&nbytes,&bit); }
    }
  s->out=b+nbytes; s->byte_off+=nbytes; s->bit_off=bit;
}

/* ========================= local urandom / autostable / gaussian ========================= */
static void
aa_urandom_i64(int64_t *vec, unsigned int nelems, int64_t mod, unsigned int log2mod,
               const uint8_t seed[32], uint64_t dom)
{
  aa_rng_state st; unsigned int k=0;
  aa_rng_init(&st, seed, dom);
  while(k < nelems){
    unsigned int l = nelems - k;
    unsigned int outlen = (log2mod * l + 7) / 8;
    uint8_t *out = (uint8_t*)malloc(outlen?outlen:1);
    unsigned int i,j;
    aa_rng_urandom(&st,out,outlen);
    for(i=0;i<l && k<nelems;i++){
      int64_t tmp=0;
      for(j=0;j<log2mod;j++){
        unsigned int qb=(i*log2mod+j)>>3;
        unsigned int rb=(i*log2mod+j)-(qb<<3);
        tmp |= (((int64_t)((out[qb]&(1u<<rb))>>rb)) << j);
      }
      if((((mod-1)-tmp) >> 63) == 0){ vec[k]=tmp; k++; }
    }
    free(out);
  }
  aa_rng_clear(&st);
}

static void
aa_urandom_mat(arith_mat_ptr M, unsigned int nrows, unsigned int ncols,
               unsigned int d, unsigned long long q, unsigned int log2q,
               const uint8_t seed[32], uint32_t dom, abdlop_arith_engine_t *e)
{
  unsigned int i,j,k; int64_t *buf=(int64_t*)malloc(d*sizeof(int64_t));
  for(i=0;i<nrows;i++)
    for(j=0;j<ncols;j++){
      uint64_t dom64 = ((uint64_t)dom << 32) | (uint64_t)(i*ncols + j);
      aa_urandom_i64(buf,d,(int64_t)q,log2q,seed,dom64);
      for(k=0;k<d;k++) e->ops->mat_set_coeff(M,i,j,k,buf[k]);
    }
  free(buf);
}

/* autostable poly: dom is passed as uint32_t widened to uint64_t (matches lazer) */
static void
aa_autostable_poly(arith_poly_ptr p, int64_t bnd, unsigned int log2,
                   const uint8_t seed[32], uint32_t dom, unsigned int d,
                   unsigned long long q, abdlop_arith_engine_t *e)
{
  unsigned int nchal=d/2, i;
  int64_t *c=(int64_t*)malloc(nchal*sizeof(int64_t));
  aa_urandom_i64(c,nchal,2*bnd+1,log2,seed,(uint64_t)dom);
  for(i=0;i<nchal;i++) e->ops->poly_set_coeff(p,i,_canonq(c[i]-bnd,q));
  e->ops->poly_set_coeff(p,nchal,0);
  for(i=nchal+1;i<d;i++) e->ops->poly_set_coeff(p,i,_canonq(-c[2*nchal-i],q));
  free(c);
}

/* gaussian sampler matching grandom.c / _grandom_sample_i32 */
typedef struct { uint64_t hi, lo; } aa_z128;
static const aa_z128 aa_CDF155[] = {
  {10894764499197476522ULL,10804844707381617341ULL},{4761708367981796450ULL,6209732027000382074ULL},
  {1476784279527800432ULL,14108379346150813303ULL},{316388870594767345ULL,17827298407763885637ULL},
  {46043503515468600ULL,18385899657892021654ULL},{4503729779335039ULL,3860889375818664979ULL},
  {294122444862326ULL,13947176349836216550ULL},{12769598070895ULL,14321894682751135119ULL},
  {367552986472ULL,10286761328368440884ULL},{7001273393ULL,2287787188970898528ULL},
  {88153536ULL,17843977990435837663ULL},{733119ULL,12174894787802692461ULL},
  {4024ULL,18067426722645776197ULL},{14ULL,10764017821655913055ULL},
  {0ULL,643125733022530080ULL},{0ULL,1014291014134832ULL},{0ULL,1055215183460ULL},
  {0ULL,724109373ULL},{0ULL,327744ULL},{0ULL,98ULL},{0ULL,0ULL} };
static inline double aa_exp_small(double x){
  double t=x*x;
  t=x-t*(1.66666666666666019037e-01+t*(-2.77777777770155933842e-03+t*(6.61375632143793436117e-05+t*(-1.65339022054652515390e-06+t*4.13813679705723846039e-08))));
  return 1.0-((x*t)/(t-2.0)-x);
}
static inline unsigned int aa_cdfsampler(aa_rng_state *st,const aa_z128 CDF[]){
  uint64_t buf[2]; unsigned int z=0; aa_rng_urandom(st,(unsigned char*)buf,16);
  while(buf[0]<=CDF[z].hi && (buf[0]<CDF[z].hi || buf[1]<CDF[z].lo)) ++z; return z;
}
static inline unsigned int aa_BerExp(aa_rng_state *st,double x){
  unsigned int b; uint64_t w[2],t,u; aa_rng_urandom(st,(unsigned char*)w,16);
  t=w[0]; u=(uint64_t)(x*(1.0/log(2.0))); x-=log(2.0)*(double)u;
  u ^= (u ^ 63ULL) & ((uint64_t)((int64_t)(63-(int64_t)u) >> 63));
  t ^= (t >> u) << u; b = 1u - (unsigned int)((t | (uint64_t)(-(int64_t)t)) >> 63);
  t=w[1]; t &= (1ULL<<53)-1; u=(uint64_t)(aa_exp_small(-x)*pow(2.0,53.0));
  b &= (unsigned int)(((int64_t)(t-u)) >> 63 ? 0 : 1);
  return b;
}
static inline int aa_gaussian155(aa_rng_state *st,double c){
  static int pos = 64; static uint64_t bits;
  int k,b; double x; const double dss=1.0/(2.0*1.55*1.55);
  do {
    if(pos>=64){ aa_rng_urandom(st,(unsigned char*)&bits,8); pos=0; }
    b=(int)(bits&1); bits>>=1; pos+=1;
    k=(int)aa_cdfsampler(st,aa_CDF155); k=(-b&(2*k))-k+b;
    x=((k-c)*(k-c)-(k-b)*(k-b))*dss;
  } while(!aa_BerExp(st,x));
  return k;
}
static void
aa_grandom_sample_i32(int32_t *ret,unsigned int nelems,unsigned int log2o,
                      const uint8_t seed[32],uint64_t dom)
{
  aa_rng_state st; unsigned int i,j; unsigned int outlen=(log2o*nelems+7)/8;
  uint8_t *out=(uint8_t*)malloc(outlen?outlen:1);
  aa_rng_init(&st,seed,dom); aa_rng_urandom(&st,out,outlen);
  for(i=0;i<nelems;i++){
    uint32_t urand=0; int32_t k; double c;
    for(j=0;j<log2o;j++){
      unsigned int qb=(i*log2o+j)>>3; unsigned int rb=(i*log2o+j)-(qb<<3);
      urand |= (((uint32_t)((out[qb]&(1u<<rb))>>rb)) << j);
    }
    c=(double)urand/(double)(1u<<log2o);
    k=aa_gaussian155(&st,c);
    ret[i]=((int32_t)k<<(int)log2o)-(int32_t)urand;
  }
  free(out); aa_rng_clear(&st);
}
/* polyvec_grandom equivalent: per element dom = (dom<<32)|(i+1) */
static void
aa_grandom_vec(arith_vec_ptr v, unsigned int nelems, unsigned int d,
               unsigned int log2o, const uint8_t seed[32], uint32_t dom,
               unsigned long long q, abdlop_arith_engine_t *e)
{
  unsigned int i,k; int32_t *samp=(int32_t*)malloc(d*sizeof(int32_t));
  for(i=0;i<nelems;i++){
    uint64_t dom64 = ((uint64_t)dom << 32) | (uint64_t)(i+1);
    aa_grandom_sample_i32(samp,d,log2o,seed,dom64);
    for(k=0;k<d;k++) e->ops->vec_set_coeff(v,i,k,_canonq((long long)samp[k],q));
  }
  free(samp);
}

/* ========================= local rejection sampling ========================= */
static void aa_mpfr_set_i128(mpfr_t f, __int128 v){
  if(v==0){ mpfr_set_ui(f,0,MPFR_RNDN); return; }
  int neg = v<0; unsigned __int128 u = neg ? (unsigned __int128)(-v) : (unsigned __int128)v;
  uint64_t hi=(uint64_t)(u>>64), lo=(uint64_t)u;
  if(hi==0){ mpfr_set_ui(f,lo,MPFR_RNDN); }
  else { mpfr_set_ui(f,hi,MPFR_RNDN); mpfr_mul_2ui(f,f,64,MPFR_RNDN); mpfr_t t; mpfr_init2(t,128); mpfr_set_ui(t,lo,MPFR_RNDN); mpfr_add(f,f,t,MPFR_RNDN); mpfr_clear(t); }
  if(neg) mpfr_neg(f,f,MPFR_RNDN);
}
static void aa_mpfr_set_limbs(mpfr_t f, const uint64_t *limbs, unsigned int nlimbs){
  mpfr_set_ui(f,0,MPFR_RNDN); mpfr_t t; mpfr_init2(t, nlimbs*64>128?(mpfr_prec_t)(nlimbs*64):128);
  int i; for(i=(int)nlimbs-1;i>=0;i--){ mpfr_mul_2ui(f,f,64,MPFR_RNDN); mpfr_set_ui(t,limbs[i],MPFR_RNDN); mpfr_add(f,f,t,MPFR_RNDN); }
  mpfr_clear(t);
}
static void aa_mul_scM_u(uint64_t out[5], const uint64_t scM[3], const uint64_t u[2]){
  unsigned __int128 acc; unsigned int i,j; for(i=0;i<5;i++) out[i]=0;
  for(i=0;i<3;i++) for(j=0;j<2;j++){
    acc=(unsigned __int128)scM[i]*u[j];
    unsigned __int128 t=(unsigned __int128)out[i+j]+(uint64_t)acc; out[i+j]=(uint64_t)t;
    unsigned __int128 c=(t>>64)+(uint64_t)(acc>>64); unsigned int k=i+j+1;
    while(c && k<5){ unsigned __int128 t2=(unsigned __int128)out[k]+(uint64_t)c; out[k]=(uint64_t)t2; c=(t2>>64)+(c>>64); k++; }
  }
}
static int
aa_rej_standard(aa_rng_state *rng, const long long *z, const long long *v,
                unsigned int nelems, const uint64_t scM[3], const uint64_t sigma2[2])
{
  uint64_t ulimbs[2]={0,0}, mu[5];
  __int128 dotzv=0, dotvv=0; unsigned int i;
  for(i=0;i<nelems;i++){ dotzv += (__int128)z[i]*( __int128)v[i]; dotvv += (__int128)v[i]*(__int128)v[i]; }
  __int128 nom128 = dotvv - 2*dotzv;
  uint64_t sigma2dbl[3];
  unsigned __int128 s2a=((unsigned __int128)sigma2[0])<<1;
  sigma2dbl[0]=(uint64_t)s2a; sigma2dbl[1]=(uint64_t)(((unsigned __int128)sigma2[1]<<1) | (s2a>>64));
  sigma2dbl[2]=(uint64_t)((unsigned __int128)sigma2[1]>>63);

  aa_rng_urandom(rng,(uint8_t*)ulimbs,16);
  aa_mul_scM_u(mu,scM,ulimbs);

  mpfr_t nom_m,denom_m,t3,t4;
  mpfr_init2(nom_m,128); mpfr_init2(denom_m,128); mpfr_init2(t3,128); mpfr_init2(t4,128);
  aa_mpfr_set_i128(nom_m,nom128);
  aa_mpfr_set_limbs(denom_m,sigma2dbl,(sigma2[1]||sigma2dbl[2])?((sigma2dbl[2])?3:2):1);
  aa_mpfr_set_limbs(t4,mu,5);
  mpfr_div(t3,nom_m,denom_m,MPFR_RNDN);
  mpfr_exp(t3,t3,MPFR_RNDN);
  mpfr_mul_2ui(t3,t3,256,MPFR_RNDN);
  int reject = (mpfr_cmp(t4,t3)>0)?1:0;
  mpfr_clears(nom_m,denom_m,t3,t4,(mpfr_ptr)0);
  return reject;
}

/* ========================= params ========================= */
void
abdlop_arith_params_set (abdlop_arith_params_t *ap,
                         unsigned long long q, unsigned int d,
                         unsigned int log2q, unsigned int D,
                         unsigned int log2m, unsigned int m1,
                         unsigned int m2, unsigned int kmsis,
                         unsigned int l, unsigned int lext,
                         long long qminus1, long long mval,
                         long long mby2v, long long gammav,
                         long long gammaby2v, long long pow2Dv,
                         long long pow2Dby2v, long long Bsqr,
                         int m_odd, unsigned int omega,
                         unsigned int log2omega,
                         unsigned int log2stdev1, int rej1,
                         const uint64_t scM1[3],
                         const uint64_t stdev1sqr[2],
                         unsigned int log2stdev2, int rej2,
                         const uint64_t scM2[3],
                         const uint64_t stdev2sqr[2])
{
  memset(ap,0,sizeof(*ap));
  ap->q=q; ap->d=d; ap->log2q=log2q; ap->D=D; ap->log2m=log2m;
  ap->log2qmd = log2q - D;
  ap->m1=m1; ap->m2=m2; ap->m22=m2-kmsis; ap->kmsis=kmsis; ap->l=l+lext;
  ap->qminus1=qminus1; ap->m=mval; ap->mby2=mby2v; ap->gamma=gammav; ap->gammaby2=gammaby2v;
  ap->pow2D=pow2Dv; ap->pow2Dby2=pow2Dby2v; ap->Bsqr=Bsqr; ap->m_odd=m_odd;
  ap->omega=omega; ap->log2omega=log2omega; ap->log2stdev1=log2stdev1; ap->log2stdev2=log2stdev2;
  ap->rej1=rej1; ap->rej2=rej2;
  memcpy(ap->scM1,scM1,sizeof(ap->scM1)); memcpy(ap->scM2,scM2,sizeof(ap->scM2));
  memcpy(ap->stdev1sqr,stdev1sqr,sizeof(ap->stdev1sqr)); memcpy(ap->stdev2sqr,stdev2sqr,sizeof(ap->stdev2sqr));
}

/* ========================= pipeline ========================= */
void
abdlop_arith_keygen (abdlop_arith_engine_t *e, const abdlop_arith_params_t *ap,
                     const uint8_t seed[32], arith_mat_ptr A1,
                     arith_mat_ptr A2prime, arith_mat_ptr Bprime)
{
  if (ap->m1 > 0) {
    aa_urandom_mat(A1, ap->kmsis, ap->m1, ap->d, ap->q, ap->log2q, seed, 0, e);
    aa_urandom_mat(A2prime, ap->kmsis, ap->m22, ap->d, ap->q, ap->log2q, seed, 1, e);
  }
  if (ap->l > 0)
    aa_urandom_mat(Bprime, ap->l, ap->m22, ap->d, ap->q, ap->log2q, seed, 2, e);
}

void
abdlop_arith_commit (abdlop_arith_engine_t *e, const abdlop_arith_params_t *ap,
                     arith_vec_ptr tA1, arith_vec_ptr tA2, arith_vec_ptr tB,
                     const arith_vec_ptr s1, const arith_vec_ptr m,
                     const arith_vec_ptr s2, const arith_mat_ptr A1,
                     const arith_mat_ptr A2prime, const arith_mat_ptr Bprime)
{
  const unsigned int m22=ap->m22, kmsis=ap->kmsis, l=ap->l, d=ap->d;
  unsigned int i,k;
  arith_vec_ptr s21 = e->ops->vec_alloc(e->ring, m22);
  for(i=0;i<m22;i++) for(k=0;k<d;k++)
    e->ops->vec_set_coeff(s21,i,k, e->ops->vec_get_coeff(s2,i,k));

  if (ap->m1 > 0) {
    arith_vec_ptr s22 = e->ops->vec_alloc(e->ring, kmsis);
    for(i=0;i<kmsis;i++) for(k=0;k<d;k++)
      e->ops->vec_set_coeff(s22,i,k, e->ops->vec_get_coeff(s2,m22+i,k));
    e->ops->vec_set(tA2, s22);
    e->ops->vec_addmul(tA2, A1, s1);
    e->ops->vec_addmul(tA2, A2prime, s21);
    for(i=0;i<kmsis;i++) for(k=0;k<d;k++){
      long long c=e->ops->vec_get_coeff(tA2,i,k), t1v,t2v;
      t1v=_p2r_c(c,&t2v,ap);
      e->ops->vec_set_coeff(tA1,i,k,t1v);
      e->ops->vec_set_coeff(tA2,i,k,t2v);
    }
    e->ops->vec_free(s22);
  }
  if (l > 0) {
    e->ops->vec_set(tB, m);
    e->ops->vec_addmul(tB, Bprime, s21);
  }
  e->ops->vec_free(s21);
}

void
abdlop_arith_hashcomm (uint8_t hash[32], abdlop_arith_engine_t *e,
                       const abdlop_arith_params_t *ap,
                       const arith_vec_ptr tA1, const arith_vec_ptr tB)
{
  const unsigned int d=ap->d, kmsis=ap->kmsis, l=ap->l;
  aa_coder cstate; aa_shake128_state hstate;
  size_t outlen = ((size_t)kmsis*d*ap->log2qmd + (size_t)l*d*ap->log2q + 7)/8 + 1;
  uint8_t *out=(uint8_t*)calloc(outlen,1); unsigned int nbits;
  aa_coder_begin(&cstate,out);
  if (ap->m1 > 0) aa_encode_urandom_vec(&cstate,e,tA1,kmsis,d,ap->log2qmd);
  if (l > 0) aa_encode_urandom_vec(&cstate,e,tB,l,d,ap->log2q);
  aa_coder_end(&cstate); nbits = aa_coder_offset(&cstate) >> 3;
  aa_shake128_init(&hstate); aa_shake128_absorb(&hstate,hash,32); aa_shake128_absorb(&hstate,out,nbits);
  aa_shake128_squeeze(&hstate,hash,32); aa_shake128_clear(&hstate);
  free(out);
}

static void
_challenge (uint8_t hash[32], abdlop_arith_engine_t *e,
            const abdlop_arith_params_t *ap, const arith_vec_ptr w1,
            arith_poly_ptr cpoly, uint8_t cseed[32])
{
  const unsigned int d=ap->d, kmsis=ap->kmsis;
  aa_coder cstate; aa_shake128_state hstate;
  size_t outlen = ((size_t)ap->log2m*d*kmsis + 7)/8 + 1;
  uint8_t *out=(uint8_t*)calloc(outlen,1); unsigned int nbits;
  aa_coder_begin(&cstate,out);
  aa_encode_urandom_vec(&cstate,e,w1,kmsis,d,ap->log2m);
  aa_coder_end(&cstate); nbits = aa_coder_offset(&cstate) >> 3;
  aa_shake128_init(&hstate); aa_shake128_absorb(&hstate,hash,32); aa_shake128_absorb(&hstate,out,nbits);
  aa_shake128_squeeze(&hstate,cseed,32); aa_shake128_clear(&hstate);
  aa_autostable_poly(cpoly,(int64_t)ap->omega,ap->log2omega,cseed,0,d,ap->q,e);
  free(out);
}

void
abdlop_arith_prove (uint8_t hash[32], arith_poly_ptr c, arith_vec_ptr z1,
                    arith_vec_ptr z21, arith_vec_ptr h,
                    const arith_vec_ptr tA2, const arith_vec_ptr s1,
                    const arith_vec_ptr s2, const arith_mat_ptr A1,
                    const arith_mat_ptr A2prime, const uint8_t seed[32],
                    abdlop_arith_engine_t *e, const abdlop_arith_params_t *ap,
                    abdlop_arith_prove_stats_t *stats)
{
  const unsigned int m1=ap->m1,m2=ap->m2,m22=ap->m22,kmsis=ap->kmsis,d=ap->d;
  aa_rng_state rngstate; uint8_t yseed[32], cseed[32]; uint32_t dom;
  unsigned int i,k,attempts=0; int rej; clock_t t_loop, t_post;

  arith_vec_ptr s21=e->ops->vec_alloc(e->ring,m22), s22=e->ops->vec_alloc(e->ring,kmsis);
  arith_vec_ptr y1v=e->ops->vec_alloc(e->ring,m1), y21v=e->ops->vec_alloc(e->ring,m22), y22v=e->ops->vec_alloc(e->ring,kmsis);
  arith_vec_ptr cs1=e->ops->vec_alloc(e->ring,m1), cs2=e->ops->vec_alloc(e->ring,m2);
  arith_vec_ptr z0m1=e->ops->vec_alloc(e->ring,m1), z0m2=e->ops->vec_alloc(e->ring,m2);
  arith_vec_ptr w=e->ops->vec_alloc(e->ring,kmsis), w1=e->ops->vec_alloc(e->ring,kmsis), w0=e->ops->vec_alloc(e->ring,kmsis), w1post=e->ops->vec_alloc(e->ring,kmsis);
  arith_poly_ptr negc=e->ops->poly_alloc(e->ring);
  long long *z1b=NULL,*cs1b=NULL,*z2b=NULL,*cs2b=NULL;
  if(ap->rej1){ z1b=(long long*)malloc(d*m1*sizeof(long long)); cs1b=(long long*)malloc(d*m1*sizeof(long long)); }
  if(ap->rej2){ z2b=(long long*)malloc(d*m2*sizeof(long long)); cs2b=(long long*)malloc(d*m2*sizeof(long long)); }

  for(i=0;i<m22;i++) for(k=0;k<d;k++) e->ops->vec_set_coeff(s21,i,k,e->ops->vec_get_coeff(s2,i,k));
  for(i=0;i<kmsis;i++) for(k=0;k<d;k++) e->ops->vec_set_coeff(s22,i,k,e->ops->vec_get_coeff(s2,m22+i,k));
  for(i=0;i<m1;i++) for(k=0;k<d;k++) e->ops->vec_set_coeff(z0m1,i,k,0);
  for(i=0;i<m2;i++) for(k=0;k<d;k++) e->ops->vec_set_coeff(z0m2,i,k,0);

  aa_rng_init(&rngstate,seed,0); aa_rng_urandom(&rngstate,yseed,32);
  dom=0; t_loop = stats?clock():0;
  while(1){
    attempts++;
    /* y1, y2 sampled fresh (overwrite) */
    aa_grandom_vec(y1v,m1,d,ap->log2stdev1,yseed,dom,ap->q,e); dom++;
    {
      /* y2 as (y21,y22) */
      arith_vec_ptr y2tmp=e->ops->vec_alloc(e->ring,m2);
      aa_grandom_vec(y2tmp,m2,d,ap->log2stdev2,yseed,dom,ap->q,e); dom++;
      for(i=0;i<m22;i++) for(k=0;k<d;k++) e->ops->vec_set_coeff(y21v,i,k,e->ops->vec_get_coeff(y2tmp,i,k));
      for(i=0;i<kmsis;i++) for(k=0;k<d;k++) e->ops->vec_set_coeff(y22v,i,k,e->ops->vec_get_coeff(y2tmp,m22+i,k));
      e->ops->vec_free(y2tmp);
    }

    /* w = y22 + A1*y1 + A2prime*y21 */
    e->ops->vec_set(w,y22v);
    e->ops->vec_addmul(w,A1,y1v);
    e->ops->vec_addmul(w,A2prime,y21v);

    /* decompose w -> w1,w0 */
    for(i=0;i<kmsis;i++) for(k=0;k<d;k++){
      long long r1,r0; _decomp_c(e->ops->vec_get_coeff(w,i,k),&r1,&r0,ap);
      e->ops->vec_set_coeff(w1,i,k,r1); e->ops->vec_set_coeff(w0,i,k,r0);
    }

    /* challenge */
    _challenge(hash,e,ap,w1,c,cseed);

    /* z1 = y1 + c*s1 ; y21 += c*s21 ; y22 += c*s22 */
    e->ops->vec_set(z1,y1v); e->ops->vec_acc_polyscalar(z1,c,s1);
    e->ops->vec_acc_polyscalar(y21v,c,s21);
    e->ops->vec_acc_polyscalar(y22v,c,s22);

    /* cs1=c*s1, cs2=c*s2 */
    e->ops->vec_set(cs1,z0m1); e->ops->vec_acc_polyscalar(cs1,c,s1);
    e->ops->vec_set(cs2,z0m2); e->ops->vec_acc_polyscalar(cs2,c,s2);

    if(ap->rej1){
      _v2bal(e,ap,z1b,z1,m1); _v2bal(e,ap,cs1b,cs1,m1);
      rej=aa_rej_standard(&rngstate,z1b,cs1b,d*m1,ap->scM1,ap->stdev1sqr);
      if(rej) continue;
    }
    if(ap->rej2){
      /* z2 = (y21||y22) after add, in balanced form */
      for(i=0;i<m22;i++) for(k=0;k<d;k++) z2b[i*d+k]=_balq(e->ops->vec_get_coeff(y21v,i,k),ap->q);
      for(i=0;i<kmsis;i++) for(k=0;k<d;k++) z2b[m22*d+i*d+k]=_balq(e->ops->vec_get_coeff(y22v,i,k),ap->q);
      _v2bal(e,ap,cs2b,cs2,m2);
      rej=aa_rej_standard(&rngstate,z2b,cs2b,d*m2,ap->scM2,ap->stdev2sqr);
      if(rej) continue;
    }

    /* y22 <- y22 - c*tA2 - w0 */
    for(k=0;k<d;k++){ long long ck=_balq(e->ops->poly_get_coeff(c,k),ap->q); e->ops->poly_set_coeff(negc,k,_canonq(-ck,ap->q)); }
    e->ops->vec_acc_polyscalar(y22v,negc,tA2);
    e->ops->vec_acc(y22v,w0,-1);

    { __int128 nrm=_l2sqr(e,y21v,m22,ap)+_l2sqr(e,y22v,kmsis,ap); if(nrm > (__int128)ap->Bsqr) continue; }
    break;
  }
  if(stats) stats->rej_ms=(double)(clock()-t_loop)/CLOCKS_PER_SEC;
  t_post = stats?clock():0;

  /* w1post = gamma*w1 - y22 */
  { arith_poly_ptr g=e->ops->poly_alloc(e->ring); e->ops->poly_set_coeff(g,0,ap->gamma);
    for(i=0;i<kmsis;i++) for(k=0;k<d;k++) e->ops->vec_set_coeff(w1post,i,k,0);
    e->ops->vec_acc_polyscalar(w1post,g,w1); e->ops->vec_acc(w1post,y22v,-1); e->ops->poly_free(g); }

  /* h = make_ghint(z=y22, r=w1post) */
  for(i=0;i<kmsis;i++) for(k=0;k<d;k++){
    long long zv=_balq(e->ops->vec_get_coeff(y22v,i,k),ap->q);
    long long rv=e->ops->vec_get_coeff(w1post,i,k);
    e->ops->vec_set_coeff(h,i,k,_make_ghint_c(rv,zv,ap));
  }

  /* z21 <- y21 */
  for(i=0;i<m22;i++) for(k=0;k<d;k++) e->ops->vec_set_coeff(z21,i,k,e->ops->vec_get_coeff(y21v,i,k));

  memcpy(hash,cseed,32);

  aa_rng_clear(&rngstate);
  e->ops->vec_free(s21); e->ops->vec_free(s22); e->ops->vec_free(y1v); e->ops->vec_free(y21v); e->ops->vec_free(y22v);
  e->ops->vec_free(cs1); e->ops->vec_free(cs2); e->ops->vec_free(z0m1); e->ops->vec_free(z0m2);
  e->ops->vec_free(w); e->ops->vec_free(w1); e->ops->vec_free(w0); e->ops->vec_free(w1post); e->ops->poly_free(negc);
  free(z1b); free(cs1b); free(z2b); free(cs2b);
  if(stats){ stats->attempts=attempts; stats->post_ms=(double)(clock()-t_post)/CLOCKS_PER_SEC; }
}

/* ========================= verify ========================= */
int
abdlop_arith_verify (uint8_t hash[32], const arith_poly_ptr c,
                     const arith_vec_ptr z1, const arith_vec_ptr z21,
                     const arith_vec_ptr h, const arith_vec_ptr tA1,
                     const arith_mat_ptr A1, const arith_mat_ptr A2prime,
                     abdlop_arith_engine_t *e, const abdlop_arith_params_t *ap)
{
  const unsigned int m1 = ap->m1;
  const unsigned int kmsis = ap->kmsis;
  const unsigned int d = ap->d;
  arith_poly_ptr negshift, c2;
  arith_vec_ptr tmp1, rem, w1v;
  uint8_t cseed[32];
  unsigned int i, k;
  int b, accept = 0;

  negshift = e->ops->poly_alloc (e->ring);
  c2 = e->ops->poly_alloc (e->ring);
  tmp1 = e->ops->vec_alloc (e->ring, kmsis);
  rem = e->ops->vec_alloc (e->ring, kmsis);
  w1v = e->ops->vec_alloc (e->ring, kmsis);

  /* -(c << D) coefficient-wise, mod q */
  for (k = 0; k < d; k++)
    {
      long long ck = _balq (e->ops->poly_get_coeff (c, k), ap->q);
      long long v = -(ck << ap->D);
      e->ops->poly_set_coeff (negshift, k, _canonq (v, ap->q));
    }

  /* tmp1 = A1*z1 + A2prime*z21 - (c<<D)*tA1 */
  e->ops->vec_addmul (tmp1, A1, z1);
  e->ops->vec_addmul (tmp1, A2prime, z21);
  e->ops->vec_acc_polyscalar (tmp1, negshift, tA1);

  /* w1 <- use_ghint(h, tmp1) */
  for (i = 0; i < kmsis; i++)
    for (k = 0; k < d; k++)
      {
        long long yv = _balq (e->ops->vec_get_coeff (h, i, k), ap->q);
        long long rv = e->ops->vec_get_coeff (tmp1, i, k);
        e->ops->vec_set_coeff (w1v, i, k, _use_ghint_c (rv, yv, ap));
      }

  /* recover challenge from w1, compare to c */
  _challenge (hash, e, ap, w1v, c2, cseed);
  b = 1;
  for (k = 0; k < d; k++)
    {
      long long a = _canonq (e->ops->poly_get_coeff (c, k), ap->q);
      long long b2 = _canonq (e->ops->poly_get_coeff (c2, k), ap->q);
      if (a != b2) { b = 0; break; }
    }
  if (!b)
    goto ret;

  /* l2(z1)^2 <= sigma1^2 * 2 * m1 * d  (2.4025*400 = 962) */
  {
    __int128 bnd = (__int128)(2 * m1 * d * 962);
    bnd <<= 2 * ap->log2stdev1;
    bnd /= 400;
    if (_l2sqr (e, z1, m1, ap) > bnd) goto ret;
  }

  /* ||tmp1 - gamma*w1||^2 <= Bsqr */
  {
    arith_poly_ptr neggamma = e->ops->poly_alloc (e->ring);
    e->ops->poly_set_coeff (neggamma, 0, _canonq (-ap->gamma, ap->q));
    e->ops->vec_set (rem, tmp1);
    e->ops->vec_acc_polyscalar (rem, neggamma, w1v);
    b = _l2sqr (e, rem, kmsis, ap) <= (__int128)ap->Bsqr;
    e->ops->poly_free (neggamma);
    if (!b) goto ret;
  }

  /* linf(h) <= m */
  if (_linf (e, h, kmsis, ap) > ap->m) goto ret;

  /* update fiat-shamir hash */
  memcpy (hash, cseed, 32);
  accept = 1;
ret:
  e->ops->poly_free (negshift);
  e->ops->poly_free (c2);
  e->ops->vec_free (tmp1);
  e->ops->vec_free (rem);
  e->ops->vec_free (w1v);
  return accept;
}
