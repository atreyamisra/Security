//see LICENSE for license
#ifndef _RISCV_RSA_ROCC_H
#define _RISCV_RSA_ROCC_H

#include "rocc.h"
#include "mmu.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>


class rsa_t : public rocc_t
{
public:
  rsa_t() {};

  const char* name() { return "rsa"; }

    typedef struct uint128 {
        uint64_t Hi;
        uint64_t Lo;
    } uint128;

  void reset()
  {
    modulusLo = 0;
    modulusHi = 0;
    privateExpLo = 0;
    privateExpHi = 0;
    pubExp = 0;
    msg_addr = 0;
    msg_len = 0;
    num=0;
    hash_addr=0;
    cLo=0;
    cHi=0;
    loAddr=0;
    hiAddr=0;
  }

  reg_t custom0(rocc_insn_t insn, reg_t xs1, reg_t xs2)
  {
    switch (insn.funct)
    {
        case 0:
          modulusLo=xs1;
          modulusHi=xs2;
          break;
        case 1:
          privateExpLo=xs1;
          privateExpHi=xs2;
          break;
        case 2:
          msg_addr=xs1;
          msg_len=xs2;
          break;
        case 3: {
          pubExp = xs1;
          char *input;
          input = (char *) malloc(msg_len * sizeof(char));
          for (uint32_t i = 0; i < msg_len; i++)
            input[i] = p->get_mmu()->load_uint8(msg_addr + i);
          printf("%s\n", input);
          num = (uint64_t) input[0]|
                (uint64_t) input[1] << 8 |
                (uint64_t) input[2] << 16 |
                (uint64_t) input[3] << 24 |
                (uint64_t) input[4] << 32 |
                (uint64_t) input[5] << 40 |
                (uint64_t) input[6] << 48 |
                (uint64_t) input[7] << 56;
          printf("NUM:%"
          PRIu64
          "\n", num);
          free(input);
          printf("%"
          PRIu64
          "----\n", modulusHi);
          printf("%"
          PRIu64
          "\n", modulusLo);
          printf("%"
          PRIu64
          "\n", privateExpLo);
          printf("%"
          PRIu64
          "\n", privateExpHi);
          printf("%"
          PRIu64
          "\n", pubExp);
          printf("%x\n", msg_addr);
          printf("%d\n", msg_len);
          uint128 modulus = {modulusHi, modulusLo};
          uint128 base = {num, 0};
          uint64_t exponent = (uint64_t) pubExp;
          uint128 c = modular_pow_reverse(base, {0,exponent}, modulus);
          cLo=c.Lo;
          cHi=c.Hi;
          printf("C%" PRIu64 "\n", c.Hi);
          printf("C%" PRIu64 "\n", c.Lo);
          break;
        }
        case 4:
          loAddr=xs1;
          hiAddr=xs2;
          p->get_mmu()->store_uint64(loAddr, cLo);
          p->get_mmu()->store_uint64(hiAddr, cHi);
          break;
        case 5:{
          loAddr=xs1;
          hiAddr=xs2;
          uint128 d = modular_pow_reverse({cHi,cLo}, {privateExpHi, privateExpLo},{modulusHi, modulusLo});
          printf("D%" PRIu64 "\n", d.Hi);
          printf("D%" PRIu64 "\n", d.Lo);
          p->get_mmu()->store_uint64(loAddr, d.Lo);
          p->get_mmu()->store_uint64(hiAddr, d.Hi);
          break;
        }

      default:
        illegal_instruction();
    }

    return -1; // accelerator currently returns nothing
  }
    static uint128 modular_pow(uint128 base, uint64_t exponent, uint128 modulus)
    {
      uint128 c = {0,1};
      for(uint64_t i=0; i<exponent; i++){
        uint128 hi = {0,0};
        uint128 lo = {0,0};
        mult128to256(base, c, hi, lo);
        //if(i%10000==0){
        //  printf("before\n");
          //printf("%" PRIu64 "\n", c.Hi);
          //printf("%" PRIu64 "\n", c.Lo);
        //}
        if((hi.Lo!=0)||(hi.Hi!=0)){
          c=mulmod(base,c,modulus);
        } else{
          uint128 quotient = {0,0};
          uint128 remainder = {0,0};
          divmod128(lo, modulus, quotient, remainder);
          c=remainder;
        }
        //if(i%10000==0){
          //printf("after\n");
          //printf("%" PRIu64 "\n", c.Hi);
          //printf("%" PRIu64 "\n", c.Lo);
        //}
      }
      return c;
    }
    static uint128 modular_pow_reverse(uint128 base, uint128 exponent, uint128 modulus)
    {
      uint128 result = {0,1};
      uint128 temp = {0,0};
      uint128 q = {0,0};
      uint128 hi = {0,0};
      uint128 lo = {0,0};
      printf("HIIIIIII%d\n", compare128(result, temp));
      printf("HIIIIIII%d\n", compare128(temp, result));
      divmod128(base, modulus, q, temp);
      base=temp;
      int i=0;
      while((exponent.Lo!=0)||(exponent.Hi!=0)){
        if((exponent.Lo%2)==1){
          mult128to256(base, result, hi, lo);
          if((hi.Lo!=0)||(hi.Hi!=0)){
            result=mulmod(base,result,modulus);
          } else{
            uint128 quotient = {0,0};
            uint128 remainder = {0,0};
            divmod128(lo, modulus, quotient, remainder);
            result=remainder;
          }
        }
        shiftright128(exponent, 1, temp);
        exponent=temp;
        temp=base;
        //printf("before\n");
        //printf("%" PRIu64 "\n", base.Hi);
        //printf("%" PRIu64 "\n", base.Lo);
        mult128to256(base, temp, hi, lo);
        if((hi.Lo!=0)||(hi.Hi!=0)){
          base=mulmod(base,temp,modulus);
        } else{
          uint128 quotient = {0,0};
          uint128 remainder = {0,0};
          divmod128(lo, modulus, quotient, remainder);
          base=remainder;
        }
        //printf("after\n");
        //printf("%" PRIu64 "\n", base.Hi);
        //printf("%" PRIu64 "\n", base.Lo);
      }
      return result;
    }
    static uint128 mulmod(uint128 a, uint128 b, uint128 modulus){
      uint128 res={0,0};
      uint128 temp={0,0};
      uint128 q={0,0};
      uint128 hi={0,0};
      uint128 lo={0,0};
      uint128 diff={0,0};
      uint128 sum={0,0};
      uint128 distance={0,0};
      temp ={18446744073709551615,18446744073709551615};
      sub128(temp, modulus, diff);
      add128(diff, {0,1}, distance);
      while((b.Lo!=0)||(b.Hi!=0)){
        if((b.Lo%2)==1){
          add128(res, a, temp);
          if(compare128(temp, res)<0){
            add128(distance, temp, sum);
            divmod128(sum, modulus, q, res);
          }else{
            divmod128(temp, modulus, q, res);
          }
        }
        mult128to256(a, {0,2}, hi, lo);
        if((hi.Lo!=0)||(hi.Hi!=0)){
          if(hi.Lo!=1)
            printf("wyf");
          divmod128(lo, modulus, q, a);
          add128(distance, a, sum);
          divmod128(sum, modulus, q, a);
        } else{
          divmod128(lo, modulus, q, a);
        }
        divmod128(b, {0,2}, q, temp);
        b=q;

      }
      return res;
    }

    static void mult128(uint128 N, uint128 M, uint128& Ans)
    {
      mult64to128(N.Lo, M.Lo, Ans.Hi, Ans.Lo);
      Ans.Hi += (N.Hi * M.Lo) + (N.Lo * M.Hi);
    }
    int static compare128(uint128 N1, uint128 N2)
    {
      return  (((N1.Hi > N2.Hi) || ((N1.Hi == N2.Hi) && (N1.Lo > N2.Lo))) ? 1 : 0)
              -  (((N1.Hi < N2.Hi) || ((N1.Hi == N2.Hi) && (N1.Lo < N2.Lo))) ? 1 : 0);
    }
    void static shiftleft128(uint128 N, unsigned S, uint128& A)
    {
      uint64_t M1, M2;
      S &= 127;

      M1 = ((((S + 127) | S) & 64) >> 6) - 1llu;
      M2 = (S >> 6) - 1llu;
      S &= 63;
      A.Hi = (N.Lo << S) & (~M2);
      A.Lo = (N.Lo << S) & M2;
      A.Hi |= ((N.Hi << S) | ((N.Lo >> (64 - S)) & M1)) & M2;

    }
    void static shiftright128(uint128 N, unsigned S, uint128& A)
    {
      uint64_t M1, M2;
      S &= 127;

      M1 = ((((S + 127) | S) & 64) >> 6) - 1llu;
      M2 = (S >> 6) - 1llu;
      S &= 63;
      A.Lo = (N.Hi >> S) & (~M2);
      A.Hi = (N.Hi >> S) & M2;
      A.Lo |= ((N.Lo >> S) | ((N.Hi << (64 - S)) & M1)) & M2;
    }
    void static inc128(uint128 N, uint128& A)
    {
      A.Lo = (N.Lo + 1);
      A.Hi = N.Hi + (((N.Lo ^ A.Lo) & N.Lo) >> 63);
    }

    void static dec128(uint128 N, uint128& A)
    {
      A.Lo = N.Lo - 1;
      A.Hi = N.Hi - (((A.Lo ^ N.Lo) & A.Lo) >> 63);
    }

    void static add128(uint128 N, uint128 M, uint128& A)
    {
      uint64_t C = (((N.Lo & M.Lo) & 1) + (N.Lo >> 1) + (M.Lo >> 1)) >> 63;
      A.Hi = N.Hi + M.Hi + C;
      A.Lo = N.Lo + M.Lo;
    }

    void static sub128(uint128 N, uint128 M, uint128& A)
    {
      A.Lo = N.Lo - M.Lo;
      uint64_t C = (((A.Lo & M.Lo) & 1) + (M.Lo >> 1) + (A.Lo >> 1)) >> 63;
      A.Hi = N.Hi - (M.Hi + C);
    }
    void static and128(uint128 N1, uint128 N2, uint128& A)
    {
      A.Hi = N1.Hi & N2.Hi;
      A.Lo = N1.Lo & N2.Lo;
    }
    size_t static ntz128(uint128 N)
    {
      return (N.Lo == 0) ? ntz64(N.Hi) + 64 : ntz64(N.Lo);
    }

    size_t static ntz64(uint64_t N)
    {
      uint64_t I = ~N;
      size_t C = ((I ^ (I + 1)) & I) >> 63;

      I = (N & 0xffffffff) + 0xffffffff;
      I = ((I & 0x100000000) ^ 0x100000000) >> 27;
      C += I;  N >>= I;

      I = (N & 0xffff) + 0xffff;
      I = ((I & 0x10000) ^ 0x10000) >> 12;
      C += I;  N >>= I;

      I = (N & 0xff) + 0xff;
      I = ((I & 0x100) ^ 0x100) >> 5;
      C += I;  N >>= I;

      I = (N & 0xf) + 0xf;
      I = ((I & 0x10) ^ 0x10) >> 2;
      C += I;  N >>= I;

      I = (N & 3) + 3;
      I = ((I & 4) ^ 4) >> 1;
      C += I;  N >>= I;

      C += ((N & 1) ^ 1);

      return C;
    }
    size_t static nlz128(uint128 N)
    {
      return (N.Hi == 0) ? nlz64(N.Lo) + 64 : nlz64(N.Hi);
    }

    size_t static nlz64(uint64_t N)
    {
      uint64_t I;
      size_t C;

      I = ~N;
      C = ((I ^ (I + 1)) & I) >> 63;

      I = (N >> 32) + 0xffffffff;
      I = ((I & 0x100000000) ^ 0x100000000) >> 27;
      C += I;  N <<= I;

      I = (N >> 48) + 0xffff;
      I = ((I & 0x10000) ^ 0x10000) >> 12;
      C += I;  N <<= I;

      I = (N >> 56) + 0xff;
      I = ((I & 0x100) ^ 0x100) >> 5;
      C += I;  N <<= I;

      I = (N >> 60) + 0xf;
      I = ((I & 0x10) ^ 0x10) >> 2;
      C += I;  N <<= I;

      I = (N >> 62) + 3;
      I = ((I & 4) ^ 4) >> 1;
      C += I;  N <<= I;

      C += (N >> 63) ^ 1;

      return C;
    }
    void static Increment(uint128& N)
    {
      uint64_t T = (N.Lo + 1);
      N.Hi += ((N.Lo ^T) & N.Lo) >> 63;
      N.Lo = T;
    }

    void static Decrement(uint128& N)
    {
      uint64_t T = (N.Lo - 1);
      N.Hi -= ((T ^ N.Lo) & T) >> 63;
      N.Lo = T;
    }

    void static Add(uint128& Ans, uint128 N, uint128 M)
    {
      uint64_t C = (((N.Lo & M.Lo) & 1) + (N.Lo >> 1) + (M.Lo >> 1)) >> 63;
      Ans.Hi = N.Hi + M.Hi + C;
      Ans.Lo = N.Lo + M.Lo;
    }

    void static Subtract(uint128& Ans, uint128 N, uint128 M)
    {
      Ans.Lo = N.Lo - M.Lo;
      uint64_t C = (((Ans.Lo & M.Lo) & 1) + (M.Lo >> 1) + (Ans.Lo >> 1)) >> 63;
      Ans.Hi = N.Hi - (M.Hi + C);
    }
    void static bindivmod128(uint128 M, uint128 N, uint128& Q, uint128 &R)
    {
      Q.Hi = Q.Lo = 0;
      size_t Shift = nlz128(N) - nlz128(M);
      shiftleft128(N, Shift, N);

      do
      {
        shiftleft128(Q, 1, Q);
        if(compare128(M, N) >= 0)
        {
          sub128(M, N, M);
          Q.Lo |= 1;
        }

        shiftright128(N, 1, N);
      }while(Shift-- != 0);

      R.Hi = M.Hi;
      R.Lo = M.Lo;
    }
    static void div128by64(const uint64_t u1, const uint64_t u0, const uint64_t v, uint64_t& q)
    {
      const uint64_t b = 1ll << 32;
      uint64_t un1, un0, vn1, vn0, q1, q0, un32, un21, un10, rhat, vs, left, right;
      size_t s;

      s = nlz64(v);
      vs = v << s;
      vn1 = vs >> 32;
      vn0 = vs & 0xffffffff;


      if (s > 0)
      {
        un32 = (u1 << s) | (u0 >> (64 - s));
        un10 = u0 << s;
      }
      else
      {
        un32 = u1;
        un10 = u0;
      }


      un1 = un10 >> 32;
      un0 = un10 & 0xffffffff;

      q1 = un32 / vn1;
      rhat = un32 % vn1;

      left = q1 * vn0;
      right = (rhat << 32) | un1;

        again1:
      if ((q1 >= b) || (left > right))
      {
        --q1;
        rhat += vn1;
        if (rhat < b)
        {
          left -= vn0;
          right = (rhat << 32) | un1;
          goto again1;
        }
      }

      un21 = (un32 << 32) + (un1 - (q1 * vs));

      q0 = un21 / vn1;
      rhat = un21 % vn1;

      left = q0 * vn0;
      right = (rhat << 32) | un0;
        again2:
      if ((q0 >= b) || (left > right))
      {
        --q0;
        rhat += vn1;
        if (rhat < b)
        {
          left -= vn0;
          right = (rhat << 32) | un0;
          goto again2;
        }
      }

      q = (q1 << 32) | q0;
    }

    void static divmod128by64(const uint64_t u1, const uint64_t u0, uint64_t v, uint64_t& q, uint64_t& r)
    {
      const uint64_t b = 1ll << 32;
      uint64_t un1, un0, vn1, vn0, q1, q0, un32, un21, un10, rhat, left, right;
      size_t s;

      s = nlz64(v);
      v <<= s;
      vn1 = v >> 32;
      vn0 = v & 0xffffffff;

      if (s > 0)
      {
        un32 = (u1 << s) | (u0 >> (64 - s));
        un10 = u0 << s;
      }
      else
      {
        un32 = u1;
        un10 = u0;
      }

      un1 = un10 >> 32;
      un0 = un10 & 0xffffffff;

      q1 = un32 / vn1;
      rhat = un32 % vn1;

      left = q1 * vn0;
      right = (rhat << 32) + un1;
        again1:
      if ((q1 >= b) || (left > right))
      {
        --q1;
        rhat += vn1;
        if (rhat < b)
        {
          left -= vn0;
          right = (rhat << 32) | un1;
          goto again1;
        }
      }

      un21 = (un32 << 32) + (un1 - (q1 * v));

      q0 = un21 / vn1;
      rhat = un21 % vn1;

      left = q0 * vn0;
      right = (rhat << 32) | un0;
        again2:
      if ((q0 >= b) || (left > right))
      {
        --q0;
        rhat += vn1;
        if (rhat < b)
        {
          left -= vn0;
          right = (rhat << 32) | un0;
          goto again2;
        }
      }

      r = ((un21 << 32) + (un0 - (q0 * v))) >> s;
      q = (q1 << 32) | q0;
    }
    static void divmod128by128(uint128 M, uint128 N, uint128& Q, uint128& R)
    {
      if (N.Hi == 0)
      {
        if (M.Hi < N.Lo)
        {
          divmod128by64(M.Hi, M.Lo, N.Lo, Q.Lo, R.Lo);
          Q.Hi = 0;
          R.Hi = 0;
          return;
        }
        else
        {
          Q.Hi = M.Hi / N.Lo;
          R.Hi = M.Hi % N.Lo;
          divmod128by64(R.Hi, M.Lo, N.Lo, Q.Lo, R.Lo);
          R.Hi = 0;
          return;
        }
      }
      else
      {
        size_t n = nlz64(N.Hi);

        uint128 v1;
        shiftleft128(N, n, v1);

        uint128 u1;
        shiftright128(M, 1, u1);

        uint128 q1;
        div128by64(u1.Hi, u1.Lo, v1.Hi, q1.Lo);
        q1.Hi = 0;
        shiftright128(q1, 63 - n, q1);

        if ((q1.Hi | q1.Lo) != 0)
        {
          dec128(q1, q1);
        }

        Q.Hi = q1.Hi;
        Q.Lo = q1.Lo;
        mult128(q1, N, q1);
        sub128(M, q1, R);

        if (compare128(R, N) >= 0)
        {
          inc128(Q, Q);
          sub128(R, N, R);
        }

        return;
      }
    }
    void static divmod128(uint128 M, uint128 N, uint128& Q, uint128& R)
    {
      size_t Nlz, Mlz, Ntz;
      int C;

      Nlz = nlz128(N);
      Mlz = nlz128(M);
      Ntz = ntz128(N);

      if(Nlz == 128)
      {
        throw 0;
      }
      else if((M.Hi | N.Hi) == 0)
      {
        Q.Hi = R.Hi = 0;
        Q.Lo = M.Lo / N.Lo;
        R.Lo = M.Lo % N.Lo;
        return;
      }
      else if(Nlz == 127)
      {
        Q = M;
        R.Hi = R.Lo = 0;
        return;
      }
      else if((Ntz + Nlz) == 127)
      {
        shiftright128(M, Ntz, Q);
        dec128(N, N);
        and128(N, M, R);
        return;
      }

      C = compare128(M, N);
      if(C < 0)
      {
        Q.Hi = Q.Lo = 0;
        R = M;
        return;
      }
      else if(C == 0)
      {
        Q.Hi = R.Hi = R.Lo = 0;
        Q.Lo = 1;
        return;
      }

      if((Nlz - Mlz) > 5)
      {
        divmod128by128(M, N, Q, R);
      }
      else
      {
        bindivmod128(M, N, Q, R);
      }
    }
    static void divmod256by128(uint128 MHi, uint128 MLo, uint128 N, uint128& R)
    {
      uint128 large={4294967296,  0};
      uint128 two={4294967296,  0};
      uint128 m1={0,  0};
      uint128 m2={0,  0};
      uint128 m3={0,  0};
      uint128 q={0,  0};
      divmod128(MHi, N, q, m1);
      divmod128(large, N, q, m2);
      divmod128(two, N, q, m3);
      printf("HereLo%" PRIu64 "\n", m1.Lo);
      printf("HereHi%" PRIu64 "\n", m1.Hi);
      printf("HereLo%" PRIu64 "\n", m2.Lo);
      printf("HereHi%" PRIu64 "\n", m2.Hi);
      printf("HereLo%" PRIu64 "\n", m3.Lo);
      printf("HereHi%" PRIu64 "\n", m3.Hi);
      uint128 hi={0,  0};
      uint128 lo={0,  0};
      mult128to256(m1, m2, hi, lo);
      printf("HereLo%" PRIu64 "\n", lo.Lo);
      printf("HereHi%" PRIu64 "\n", lo.Hi);
      if((hi.Lo!=0)||(hi.Hi!=0)){
        uint128 lo1 = {0,0};
        divmod256by128(hi, lo, N, lo1);
        lo=lo1;
      }
      printf("HereLo%" PRIu64 "\n", lo.Lo);
      printf("HereHi%" PRIu64 "\n", lo.Hi);
      mult128to256(m3, lo, hi, lo);
      if((hi.Lo!=0)||(hi.Hi!=0)){
        uint128 lo1 = {0,0};
        divmod256by128(hi, lo, N, lo1);
        lo=lo1;
      }
      printf("HereLo%" PRIu64 "\n", lo.Lo);
      printf("HereHi%" PRIu64 "\n", lo.Hi);
      uint128 r1={0,  0};
      uint128 r2={0,  0};
      divmod128(lo,N,q,r1);
      uint128 result={0,  0};
      divmod128(MLo, N, q, r2);
      add128(r1,r2,result);
      printf("HereLo%" PRIu64 "\n", r1.Lo);
      printf("HereHi%" PRIu64 "\n", r1.Hi);
      printf("HereLo%" PRIu64 "\n", r2.Lo);
      printf("HereHi%" PRIu64 "\n", r2.Hi);
      divmod128(result, N, q, R);
    }

    //Now for multiply
    void static mult64to128(uint64_t u, uint64_t v, uint64_t& h, uint64_t& l)
    {
      uint64_t u1 = (u & 0xffffffff);
      uint64_t v1 = (v & 0xffffffff);
      uint64_t t = (u1 * v1);
      uint64_t w3 = (t & 0xffffffff);
      uint64_t k = (t >> 32);

      u >>= 32;
      t = (u * v1) + k;
      k = (t & 0xffffffff);
      uint64_t w1 = (t >> 32);

      v >>= 32;
      t = (u1 * v) + k;
      k = (t >> 32);

      h = (u * v) + w1 + k;
      l = (t << 32) + w3;
    }
    void static mult128to256(uint128 N, uint128 M, uint128& H, uint128& L)
    {
      mult64to128(N.Hi, M.Hi, H.Hi, H.Lo);
      mult64to128(N.Lo, M.Lo, L.Hi, L.Lo);

      uint128 T;
      mult64to128(N.Hi, M.Lo, T.Hi, T.Lo);
      L.Hi += T.Lo;
      if(L.Hi < T.Lo)  // if L.Hi overflowed
      {
        Increment(H);
      }
      H.Lo += T.Hi;
      if(H.Lo < T.Hi)  // if H.Lo overflowed
      {
        ++H.Hi;
      }

      mult64to128(N.Lo, M.Hi, T.Hi, T.Lo);
      L.Hi += T.Lo;
      if(L.Hi < T.Lo)  // if L.Hi overflowed
      {
        Increment(H);
      }
      H.Lo += T.Hi;
      if(H.Lo < T.Hi)  // if H.Lo overflowed
      {
        ++H.Hi;
      }
    }
private:
    reg_t modulusLo;
    reg_t modulusHi;
    reg_t privateExpLo;
    reg_t privateExpHi;
    reg_t pubExp;
    reg_t msg_addr;
    reg_t msg_len;
    uint64_t num;
    reg_t hash_addr;
    uint64_t cLo;
    uint64_t cHi;
    reg_t loAddr;
    reg_t hiAddr;
};
REGISTER_EXTENSION(rsa, []() { return new rsa_t; })
#endif
