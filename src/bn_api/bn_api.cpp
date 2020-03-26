#include "bn_api.h"

#include <iostream>
#include <iomanip>

using namespace wabt;
using namespace wabt::interp;

typedef unsigned __int128 uint128_t;

const bool tracing = false;

wabt::interp::Memory* GlobalMem = nullptr;

intx::uint256 BignumOne = intx::from_string<intx::uint256>("1");

intx::uint256 two_pow_256 = intx::from_string<intx::uint256>("115792089237316195423570985008687907853269984665640564039457584007913129639935");

// Mask128 = 0xffffffffffffffffffffffffffffffff
intx::uint128 Mask128 = intx::from_string<intx::uint128>("340282366920938463463374607431768211455");


intx::uint256 FqModulus = intx::from_string<intx::uint256>("21888242871839275222246405745257275088696311157297823662689037894645226208583");
intx::uint256 FqInv = intx::from_string<intx::uint256>("211173256549385567650468519415768310665");
intx::uint256 FqRsquared = intx::from_string<intx::uint256>("3096616502983703923843567936837374451735540968419076528771170197431451843209");



/*** Fq constants here are for secp256k1 ***/
// modulus = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f
// inv = 0xbcb223fedc24a059d838091dd2253531
// r_squared = 0x1000007a2000e90a1
//intx::uint256 FqModulus = intx::from_string<intx::uint256>("115792089237316195423570985008687907853269984665640564039457584007908834671663");
//intx::uint256 FqInv = intx::from_string<intx::uint256>("250819822124803770581580479000962479409");
//intx::uint256 FqRsquared = intx::from_string<intx::uint256>("18446752466076602529");


/** Fr constants are for bn128. leaving them here just so the code compiles **/
//intx::uint256 FrModulus = intx::from_string<intx::uint256>("21888242871839275222246405745257275088548364400416034343698204186575808495617");
//intx::uint256 FrInv = intx::from_string<intx::uint256>("134950519161967129512891662437158223871");
//intx::uint256 FrRsquared = intx::from_string<intx::uint256>("944936681149208446651664254269745548490766851729442924617792859073125903783");


void trace_words(uint8_t *mem, size_t len) {
    std::cout << std::hex << std::setw(2) << std::setfill('0');
    for (auto i = 0; i < len; i++) {
        for (auto j = 0; j < 32; j++) {
            std::cout << static_cast<int>(*(mem + j));
        }
        std::cout << std::dec << std::endl << std::hex;
    }

    std::cout << std::dec << std::endl;
}

static void trace_word(uint8_t *mem) {
    std::cout << std::hex << std::setw(2) << std::setfill('0');
    for (auto j = 0; j < 32; j++) {
        std::cout << static_cast<int>(*(mem + j));
    }

    std::cout << std::dec << std::endl;
}


// non-interleaved (slow)
void montgomery_multiplication_256(uint64_t* x, uint64_t* y, uint64_t* m, uint64_t* inv_ptr, uint64_t* out) {
  using intx::uint512;

  intx::uint256* a = reinterpret_cast<intx::uint256*>(x);
  intx::uint256* b = reinterpret_cast<intx::uint256*>(y);
  intx::uint256* mod = reinterpret_cast<intx::uint256*>(m);
  intx::uint256* inv = reinterpret_cast<intx::uint256*>(inv_ptr);
  intx::uint256* ret = reinterpret_cast<intx::uint256*>(out);

  //std::cout << "montgomery_multiplication_256 using non-interleaved.  a: " << intx::to_string(*a) << "  b: " << intx::to_string(*b) << std::endl;

  auto res1 = uint512{*a} * uint512{*b};
  //auto k0 = ((inv * res1).lo).lo;
  auto k0 = (uint512{*inv} * res1).lo & Mask128;
  auto res2 = ((uint512{k0} * uint512{*mod}) + res1) >> 128;
  auto k1 = (res2 * uint512{*inv}).lo & Mask128;
  auto result = ((uint512{k1} * uint512{*mod}) + res2) >> 128;
  if (result >= *mod) {
    result = result - *mod;
  }

  intx::uint256 result_256 = result.lo;
  //std::cout << "montgomery_multiplication_256 using non-interleaved.  result: " << intx::to_string(result_256) << std::endl;
  intx::uint128 result_low_128 = result_256.lo;
  intx::uint128 result_high_128 = result_256.hi;

  out[0] = result_low_128.lo;
  out[1] = result_low_128.hi;
  out[2] = result_high_128.lo;
  out[3] = result_high_128.hi;
}



// interleaved (fast)
void montgomery_multiplication_256_interleaved(uint64_t* x, uint64_t* y, uint64_t* m, uint64_t* inv, uint64_t *out){

  uint64_t A[] = {0,0,0,0,0,0,0,0,0}; // we need 9 64-bit limbs, the 9th limb in case x*y (before subtracting the modulus) is greater than 256 bits
  for (int i=0; i<4; i++){
    uint64_t ui = (A[i]+x[i]*y[0])*inv[0];
    uint64_t carry = 0;
#pragma unroll
    for (int j=0; j<4; j++){
      uint128_t xiyj = (uint128_t)x[i]*y[j];
      uint128_t uimj = (uint128_t)ui*m[j];
      uint128_t partial_sum = xiyj+carry;
      uint128_t sum = uimj+A[i+j]+partial_sum;
      A[i+j] = (uint64_t)sum;
      carry = sum>>64;

      if (sum<partial_sum){
        int k=2;
        while (i+j+k<8 && A[i+j+k]==(uint64_t)0-1 ){
          A[i+j+k]=0;
          k++;
        }
        if (i+j+k<9)
          A[i+j+k]+=1;
      }

    }
    A[i+4]+=carry;
  }

  // copy A[4:7] to out
  for (int i=0; i < 256/64; i++)
    out[i] = A[i + 256/64];

  uint64_t outMulOver256 = A[8];

  int geq = 1; // out >= modulus
  if (outMulOver256 > 0) {
    // out is greater than 256 bits
    geq = 1;
  } else {
    // out is 256 bits or less, so compare to modulus
    for (int i=256/64 - 1;i>=0;i--){
      if (out[i]<m[i]){
        // if outMul[i] is less than m[i], then out is less than the modulus and we don't need to subtract.
        geq = 0;
        break;
      } else if (out[i]>m[i]) {
        // if outMul is greater than m[i], then we need to subtract the modulus from outMul
        geq = 1;
        break;
      }
      // if this limb is equal for m[i] and outMul[i], go to next limb [i--]
    }
  }

  //if out >= modulus, then out = out - modulus
  if (geq){
    uint64_t carry=0;
#pragma unroll
    for (int i=0; i<4;i++){

      // if a limb subtraction gets a negative, then set the carry bit to 1
      // 0x0000000000000f000 - 0x0000000000a00000 = -10424320 = 0xffffffffff60f000

      // uint64_t mod_temp = m[i]-carry; // this doesn't work
      // if m is 0x..ffff and carry is 1, then mod_temp is 0x..fffe, which gets interpreted as [outMul[i] - -2] == [outMul[i] + 2 when subrtracted like [outMul[i] - mod_temp]

      // this way works
      uint64_t out_temp = out[i] - m[i] - carry;

      // if out[i] is larger than m[i], then carry should be 0
      // if out[i] is less than m[i], then subtracting m[i] will result in a carry.
      // if m[i] is less than carry, then subtracting m[i]-carry will result in another carry
      carry = (out[i]<m[i] || m[i]<carry) ? 1:0;

      // write result to return memory offset
      out[i] = out_temp;
    }
    // on i=4, if out[4] is 1, then subtracting the modulus will lead to subtracting the last carry bit from out[4]
  }
}

/*

// interleaved (broken)

void montgomery_multiplication_256(uint64_t* const x, uint64_t* const y, uint64_t* const m, uint64_t* const inv, uint64_t *out){
  uint64_t A[256/64*2] = {0};
  for (int i=0; i<256/64; i++){
    uint64_t ui = (A[i]+x[i]*y[0])*inv[0];
    uint64_t carry = 0;
#pragma unroll
    for (int j=0; j<256/64; j++){
      __uint128_t xiyj = (__uint128_t)x[i]*y[j];
      __uint128_t uimj = (__uint128_t)ui*m[j];
      __uint128_t partial_sum = xiyj+carry+A[i+j];
      __uint128_t sum = uimj+partial_sum;
      A[i+j] = (uint64_t)sum;
      carry = sum>>64;

      if (sum<partial_sum){
        int k=2;
        while ( i+j+k<256/64*2 && A[i+j+k]==(uint64_t)0-1 ){
          A[i+j+k]=0;
          k++;
        }
        if (i+j+k<9)
          A[i+j+k]+=1;
      }

    }
    A[i+256/64]+=carry;
  }
  for (int i=0; i<=256/64;i++)
    out[i] = A[i+256/64];
  // check if m <= out
  int leq = 1;
  for (int i=256/64;i>=0;i--){
    if (m[i]>out[i]){
      leq = 0;
      break;
    }
    else if (m[i]<out[i]){
      break;
    }
  }
  // if leq, then perform final subtraction
  if (leq){
    uint64_t carry=0;
#pragma unroll
    for (int i=0; i<=256/64;i++){
      uint64_t temp = m[i]-carry;
      out[i] = temp-out[i];
      carry = (temp<out[i] || m[i]<carry) ? 1:0;
    }
  }

}
*/

/*
void BNAPI::SetMemory(wabt::interp::Memory *memory) {
    this->memory = memory;
}
*/

void BNAPI::AddHostFunctions(wabt::interp::Environment* env, wabt::interp::HostModule *host_module) {


  //wabt::interp::Memory* mem = env->GetMemory(0);
  //this->memory = mem;
  GlobalMem = env->GetMemory(0);
  //wabt::interp::Memory* this->memory = env->GetMemory(0);



	host_module->AppendFuncExport("mulmodmont256", {{Type::I32, Type::I32, Type::I32, Type::I32, Type::I32}, {}}, 
        [env]( const interp::HostFunc* host_func, 
             const interp::FuncSignature *signature,
             const interp::TypedValues &args,
             interp::TypedValues &results) {

            printf("bn_api mulmodmont256.\n");

            wabt::interp::Memory* mem = env->GetMemory(0);

            uint32_t arg1_offset = static_cast<uint32_t>(args[0].value.i32);
            uint32_t arg2_offset = static_cast<uint32_t>(args[1].value.i32);
            uint32_t mod_offset = static_cast<uint32_t>(args[2].value.i32);
            uint32_t inv_offset = static_cast<uint32_t>(args[3].value.i32);
            uint32_t ret_offset = static_cast<uint32_t>(args[4].value.i32);


            //this->mulmodmont256(arg1_offset, arg2_offset, mod_offset, inv_offset, ret_offset);
            //return interp::ResultType::Ok;
            //printf("bn_api mulmodmont256 returning.\n");

            
            uint64_t* a = reinterpret_cast<uint64_t*>(&(mem->data[arg1_offset]));
            uint64_t* b = reinterpret_cast<uint64_t*>(&(mem->data[arg2_offset]));
            uint64_t* modulus = reinterpret_cast<uint64_t*>(&(mem->data[mod_offset]));
            uint64_t* inv = reinterpret_cast<uint64_t*>(&(mem->data[inv_offset]));
            uint64_t* ret = reinterpret_cast<uint64_t*>(&(mem->data[ret_offset]));

            montgomery_multiplication_256(a, b, modulus, inv, ret);
            return interp::Result::Ok;
        });


	host_module->AppendFuncExport("addmod256", {{Type::I32, Type::I32, Type::I32, Type::I32}, {}},
        [env]( const interp::HostFunc* host_func, 
             const interp::FuncSignature *signature,
             const interp::TypedValues &args,
             interp::TypedValues &results) {
            
            //printf("bn_api addmod256.\n");

            wabt::interp::Memory* mem = env->GetMemory(0);

            uint32_t arg1_offset = static_cast<uint32_t>(args[0].value.i32);
            uint32_t arg2_offset = static_cast<uint32_t>(args[1].value.i32);
            uint32_t mod_offset = static_cast<uint32_t>(args[2].value.i32);
            uint32_t ret_offset = static_cast<uint32_t>(args[3].value.i32);

            intx::uint256* a = reinterpret_cast<intx::uint256*>(&(mem->data[arg1_offset]));
          	intx::uint256* b = reinterpret_cast<intx::uint256*>(&(mem->data[arg2_offset]));
            intx::uint256* modulus = reinterpret_cast<intx::uint256*>(&(mem->data[mod_offset]));
          	intx::uint256* ret_mem = reinterpret_cast<intx::uint256*>(&(mem->data[ret_offset]));

            intx::uint512 ret_full = intx::uint512{0,*a} + intx::uint512{0,*b};

          	if (ret_full >= *modulus) {
          	  ret_full -= intx::uint512{0, *modulus};
          	}

          	*ret_mem = ret_full.lo;

            //this->addmod256(arg1_offset, arg2_offset, mod_offset, ret_offset);

            //printf("bn_api addmod256 returning.\n");

            //return interp::ResultType::Ok;
            return interp::Result::Ok;
        });

	host_module->AppendFuncExport("submod256", {{Type::I32, Type::I32, Type::I32, Type::I32}, {}},
        [env]( const interp::HostFunc* host_func, 
             const interp::FuncSignature *signature,
             const interp::TypedValues &args,
             interp::TypedValues &results) {

            //printf("bn_api submod256.\n");

            uint32_t arg1_offset = static_cast<uint32_t>(args[0].value.i32);
            uint32_t arg2_offset = static_cast<uint32_t>(args[1].value.i32);
            uint32_t mod_offset = static_cast<uint32_t>(args[2].value.i32);
            uint32_t ret_offset = static_cast<uint32_t>(args[3].value.i32);


            /* ** using this->memory
            intx::uint256* a = reinterpret_cast<intx::uint256*>(&(this->memory->data[arg1_offset]));
            intx::uint256* b = reinterpret_cast<intx::uint256*>(&(this->memory->data[arg2_offset]));
            intx::uint256* modulus = reinterpret_cast<intx::uint256*>(&(this->memory->data[mod_offset]));
            intx::uint256* ret_mem = reinterpret_cast<intx::uint256*>(&(this->memory->data[ret_offset]));
            
            this->submod256(arg1_offset, arg2_offset, mod_offset, ret_offset);
            */
            

            wabt::interp::Memory* mem = env->GetMemory(0);

            intx::uint256* a = reinterpret_cast<intx::uint256*>(&(mem->data[arg1_offset]));
            intx::uint256* b = reinterpret_cast<intx::uint256*>(&(mem->data[arg2_offset]));
            intx::uint256* modulus = reinterpret_cast<intx::uint256*>(&(mem->data[mod_offset]));
            intx::uint256* ret_mem = reinterpret_cast<intx::uint256*>(&(mem->data[ret_offset]));

            if (*a < *b) {
              *ret_mem = (*a + *modulus) - *b;
            } else {
              *ret_mem = *a - *b;
            }


            /*
            if (tracing) {
                std::cout << " = " << intx::to_string(*ret_mem) << std::endl;
            }
            */

            //printf("bn_api submod256 returning.\n");

            //return interp::ResultType::Ok;
            return interp::Result::Ok;
        });




    host_module->AppendFuncExport("bignum_f1m_mul", {{Type::I32, Type::I32, Type::I32}, {}}, 
          [env]( const interp::HostFunc* host_func, 
               const interp::FuncSignature *signature,
               const interp::TypedValues &args,
               interp::TypedValues &results) {

              uint32_t arg1_offset = static_cast<uint32_t>(args[0].value.i32);
              uint32_t arg2_offset = static_cast<uint32_t>(args[1].value.i32);
              uint32_t ret_offset = static_cast<uint32_t>(args[2].value.i32);

              /*
              wabt::interp::Memory* mem = env->GetMemory(0);
              intx::uint256* a = reinterpret_cast<intx::uint256*>(&(mem->data[arg1_offset]));
              intx::uint256* b = reinterpret_cast<intx::uint256*>(&(mem->data[arg2_offset]));
              //intx::uint256* modulus = reinterpret_cast<intx::uint256*>(&(mem->data[mod_offset]));
              intx::uint256* ret_mem = reinterpret_cast<intx::uint256*>(&(mem->data[ret_offset]));
              */

              /*
              uint64_t* a = reinterpret_cast<uint64_t*>(&(this->memory->data[a_offset]));
              uint64_t* b = reinterpret_cast<uint64_t*>(&(this->memory->data[b_offset]));
              uint64_t* ret = reinterpret_cast<uint64_t*>(&(this->memory->data[ret_offset]));
              */

              wabt::interp::Memory* mem = env->GetMemory(0);
              uint64_t* a = reinterpret_cast<uint64_t*>(&(mem->data[arg1_offset]));
              uint64_t* b = reinterpret_cast<uint64_t*>(&(mem->data[arg2_offset]));
              uint64_t* ret = reinterpret_cast<uint64_t*>(&(mem->data[ret_offset]));


              //mul_mod(a, b, &FqModulus, &FqInv, ret);
              //montgomery_multiplication_256(a, b, &FqModulus, &FqInv, ret);


              uint64_t *mod = reinterpret_cast<uint64_t *>(&(FqModulus));
              uint64_t *inv = reinterpret_cast<uint64_t *>(&(FqInv));

              montgomery_multiplication_256(a, b, mod, inv, ret);


              //this->f1m_mul(arg1_offset, arg2_offset, ret_offset);
              return interp::Result::Ok;
          });


          host_module->AppendFuncExport("bignum_f1m_mul", {{Type::I32, Type::I32, Type::I32}, {}}, 
                [env]( const interp::HostFunc* host_func, 
                     const interp::FuncSignature *signature,
                     const interp::TypedValues &args,
                     interp::TypedValues &results) {

                    uint32_t arg1_offset = static_cast<uint32_t>(args[0].value.i32);
                    uint32_t arg2_offset = static_cast<uint32_t>(args[1].value.i32);
                    uint32_t ret_offset = static_cast<uint32_t>(args[2].value.i32);

                    wabt::interp::Memory* mem = env->GetMemory(0);
                    uint64_t* a = reinterpret_cast<uint64_t*>(&(mem->data[arg1_offset]));
                    uint64_t* b = reinterpret_cast<uint64_t*>(&(mem->data[arg2_offset]));
                    uint64_t* ret = reinterpret_cast<uint64_t*>(&(mem->data[ret_offset]));

                    uint64_t *mod = reinterpret_cast<uint64_t *>(&(FqModulus));
                    uint64_t *inv = reinterpret_cast<uint64_t *>(&(FqInv));

                    montgomery_multiplication_256(a, b, mod, inv, ret);

                    return interp::Result::Ok;
                });


    host_module->AppendFuncExport("bignum_f1m_square", {{Type::I32, Type::I32}, {}},
          [env]( const interp::HostFunc* host_func, 
               const interp::FuncSignature *signature,
               const interp::TypedValues &args,
               interp::TypedValues &results) {

               uint32_t arg1_offset = static_cast<uint32_t>(args[0].value.i32);
               uint32_t ret_offset = static_cast<uint32_t>(args[1].value.i32);

               wabt::interp::Memory* mem = env->GetMemory(0);
               uint64_t* a = reinterpret_cast<uint64_t*>(&(mem->data[arg1_offset]));
               uint64_t* ret = reinterpret_cast<uint64_t*>(&(mem->data[ret_offset]));

               uint64_t *mod = reinterpret_cast<uint64_t *>(&(FqModulus));
               uint64_t *inv = reinterpret_cast<uint64_t *>(&(FqInv));

               montgomery_multiplication_256(a, a, mod, inv, ret);

               return interp::Result::Ok;
          });


    host_module->AppendFuncExport("bignum_f1m_add", {{Type::I32, Type::I32, Type::I32}, {}},
          [env]( const interp::HostFunc* host_func, 
               const interp::FuncSignature *signature,
               const interp::TypedValues &args,
               interp::TypedValues &results) {

               uint32_t arg1_offset = static_cast<uint32_t>(args[0].value.i32);
               uint32_t arg2_offset = static_cast<uint32_t>(args[1].value.i32);
               uint32_t ret_offset = static_cast<uint32_t>(args[2].value.i32);

               wabt::interp::Memory* mem = env->GetMemory(0);
               intx::uint256* a = reinterpret_cast<intx::uint256*>(&(mem->data[arg1_offset]));
               intx::uint256* b = reinterpret_cast<intx::uint256*>(&(mem->data[arg2_offset]));
               intx::uint256* ret_mem = reinterpret_cast<intx::uint256*>(&(mem->data[ret_offset]));

               intx::uint512 ret_full = intx::uint512{0,*a} + intx::uint512{0,*b};

               if (ret_full >= FqModulus) {
                 ret_full -= intx::uint512{0, FqModulus};
               }

               *ret_mem = ret_full.lo;


               return interp::Result::Ok;

          });

    host_module->AppendFuncExport("bignum_f1m_sub", {{Type::I32, Type::I32, Type::I32}, {}},
          [env]( const interp::HostFunc* host_func, 
               const interp::FuncSignature *signature,
               const interp::TypedValues &args,
               interp::TypedValues &results) {

               uint32_t arg1_offset = static_cast<uint32_t>(args[0].value.i32);
               uint32_t arg2_offset = static_cast<uint32_t>(args[1].value.i32);
               uint32_t ret_offset = static_cast<uint32_t>(args[2].value.i32);

               wabt::interp::Memory* mem = env->GetMemory(0);
               intx::uint256* a = reinterpret_cast<intx::uint256*>(&(mem->data[arg1_offset]));
               intx::uint256* b = reinterpret_cast<intx::uint256*>(&(mem->data[arg2_offset]));
               intx::uint256* ret_mem = reinterpret_cast<intx::uint256*>(&(mem->data[ret_offset]));

               if (*a < *b) {
                  *ret_mem = (*a + FqModulus) - *b;
                } else {
                  *ret_mem = *a - *b;
                }

                return interp::Result::Ok;

          });


    host_module->AppendFuncExport("bignum_f1m_fromMontgomery", {{Type::I32, Type::I32}, {}},
          [env]( const interp::HostFunc* host_func, 
               const interp::FuncSignature *signature,
               const interp::TypedValues &args,
               interp::TypedValues &results) {

                 uint32_t arg1_offset = static_cast<uint32_t>(args[0].value.i32);
                 uint32_t ret_offset = static_cast<uint32_t>(args[1].value.i32);

                 wabt::interp::Memory* mem = env->GetMemory(0);
                 uint64_t* a = reinterpret_cast<uint64_t*>(&(mem->data[arg1_offset]));
                 uint64_t* ret = reinterpret_cast<uint64_t*>(&(mem->data[ret_offset]));

                 uint64_t *mod = reinterpret_cast<uint64_t *>(&(FqModulus));
                 uint64_t *inv = reinterpret_cast<uint64_t *>(&(FqInv));


                 //fromMont(a, &FqModulus, &FqInv, ret);
                 uint64_t* b = reinterpret_cast<uint64_t*>(&BignumOne);
                 //uint64_t *mod = reinterpret_cast<uint64_t *>(modulus);
                 //uint64_t *inv = reinterpret_cast<uint64_t *>(inverse);

                 montgomery_multiplication_256(a, b, mod, inv, ret);

                 return interp::Result::Ok;
          });

    host_module->AppendFuncExport("bignum_f1m_toMontgomery", {{Type::I32, Type::I32}, {}},
          [env]( const interp::HostFunc* host_func, 
               const interp::FuncSignature *signature,
               const interp::TypedValues &args,
               interp::TypedValues &results) {

              uint32_t arg1_offset = static_cast<uint32_t>(args[0].value.i32);
              uint32_t ret_offset = static_cast<uint32_t>(args[1].value.i32);

              wabt::interp::Memory* mem = env->GetMemory(0);
              uint64_t* a = reinterpret_cast<uint64_t*>(&(mem->data[arg1_offset]));
              uint64_t* ret = reinterpret_cast<uint64_t*>(&(mem->data[ret_offset]));

              uint64_t *mod = reinterpret_cast<uint64_t *>(&(FqModulus));
              uint64_t *inv = reinterpret_cast<uint64_t *>(&(FqInv));
              uint64_t *r_squared = reinterpret_cast<uint64_t *>(&(FqRsquared));

              //fromMont(a, &FqModulus, &FqInv, ret);
              //uint64_t* b = reinterpret_cast<uint64_t*>(&BignumOne);

              //uint64_t *mod = reinterpret_cast<uint64_t *>(modulus);
              //uint64_t *inv = reinterpret_cast<uint64_t *>(inverse);

              montgomery_multiplication_256(a, r_squared, mod, inv, ret);

              return interp::Result::Ok;
          });



  	host_module->AppendFuncExport("bignum_int_mul", {{Type::I32, Type::I32, Type::I32}, {}},
          [env]( const interp::HostFunc* host_func, 
               const interp::FuncSignature *signature,
               const interp::TypedValues &args,
               interp::TypedValues &results) {

                 uint32_t arg1_offset = static_cast<uint32_t>(args[0].value.i32);
                 uint32_t arg2_offset = static_cast<uint32_t>(args[1].value.i32);
                 uint32_t ret_offset = static_cast<uint32_t>(args[2].value.i32);

                 wabt::interp::Memory* mem = env->GetMemory(0);
                 intx::uint256* a = reinterpret_cast<intx::uint256*>(&(mem->data[arg1_offset]));
                 intx::uint256* b = reinterpret_cast<intx::uint256*>(&(mem->data[arg2_offset]));
                 intx::uint256* ret_mem = reinterpret_cast<intx::uint256*>(&(mem->data[ret_offset]));

                 *ret_mem = (*a * *b ) % two_pow_256;

                  return interp::Result::Ok;
          });

  	host_module->AppendFuncExport("bignum_int_add", {{Type::I32, Type::I32, Type::I32}, {Type::I32}},
          [env]( const interp::HostFunc* host_func, 
               const interp::FuncSignature *signature,
               const interp::TypedValues &args,
               interp::TypedValues &results) {

              uint32_t arg1_offset = static_cast<uint32_t>(args[0].value.i32);
              uint32_t arg2_offset = static_cast<uint32_t>(args[1].value.i32);
              uint32_t ret_offset = static_cast<uint32_t>(args[2].value.i32);

              wabt::interp::Memory* mem = env->GetMemory(0);
              intx::uint256* a = reinterpret_cast<intx::uint256*>(&(mem->data[arg1_offset]));
              intx::uint256* b = reinterpret_cast<intx::uint256*>(&(mem->data[arg2_offset]));
              intx::uint256* ret_mem = reinterpret_cast<intx::uint256*>(&(mem->data[ret_offset]));

              const auto add_res = add_with_carry(*a, *b);

              *ret_mem = std::get<0>(add_res);

              uint32_t carry = 0;
              if (std::get<1>(add_res) > 0) {
                carry = 1;
              }

              //return carry;
              results[0].set_i32(carry);

              //uint32_t result = this->add256(arg1_offset, arg2_offset, ret_offset);

              return interp::Result::Ok;
          });

  	host_module->AppendFuncExport("bignum_int_sub", {{Type::I32, Type::I32, Type::I32}, {Type::I32}},
          [env]( const interp::HostFunc* host_func, 
               const interp::FuncSignature *signature,
               const interp::TypedValues &args,
               interp::TypedValues &results) {

                 uint32_t arg1_offset = static_cast<uint32_t>(args[0].value.i32);
                 uint32_t arg2_offset = static_cast<uint32_t>(args[1].value.i32);
                 uint32_t ret_offset = static_cast<uint32_t>(args[2].value.i32);

                 wabt::interp::Memory* mem = env->GetMemory(0);
                 intx::uint256* a = reinterpret_cast<intx::uint256*>(&(mem->data[arg1_offset]));
                 intx::uint256* b = reinterpret_cast<intx::uint256*>(&(mem->data[arg2_offset]));
                 intx::uint256* ret_mem = reinterpret_cast<intx::uint256*>(&(mem->data[ret_offset]));

                 uint32_t carry = 0;

               	if (*a < *b) {
               	  carry = 1;
               	}

               	*ret_mem = *a - *b;

                 //return carry;
                 results[0].set_i32(carry);

                 return interp::Result::Ok;
          });

  	host_module->AppendFuncExport("bignum_int_div", {{Type::I32, Type::I32, Type::I32, Type::I32}, {}},
          [env]( const interp::HostFunc* host_func, 
               const interp::FuncSignature *signature,
               const interp::TypedValues &args,
               interp::TypedValues &results) {

              uint32_t arg1_offset = static_cast<uint32_t>(args[0].value.i32);
              uint32_t arg2_offset = static_cast<uint32_t>(args[1].value.i32);
              uint32_t carry_offset = static_cast<uint32_t>(args[2].value.i32);
              uint32_t result_offset = static_cast<uint32_t>(args[3].value.i32);

              wabt::interp::Memory* mem = env->GetMemory(0);
              intx::uint256* a = reinterpret_cast<intx::uint256*>(&(mem->data[arg1_offset]));
              intx::uint256* b = reinterpret_cast<intx::uint256*>(&(mem->data[arg2_offset]));
              intx::uint256* ret_remainder_mem = reinterpret_cast<intx::uint256*>(&(mem->data[result_offset]));
              intx::uint256* ret_quotient_mem = reinterpret_cast<intx::uint256*>(&(mem->data[carry_offset]));

              const auto div_res = udivrem(*a, *b);
              *ret_quotient_mem = div_res.quot;
              *ret_remainder_mem = div_res.rem;

              //this->div256(arg1_offset, arg2_offset, carry_offset, result_offset);

              return interp::Result::Ok;
          });


}



/*
void BNAPI::mul256(uint32_t &a_offset, uint32_t &b_offset, uint32_t &ret_offset) {
	intx::uint256* a = reinterpret_cast<intx::uint256*>(&(this->memory->data[a_offset]));
	intx::uint256* b = reinterpret_cast<intx::uint256*>(&(this->memory->data[b_offset]));
	intx::uint256* ret_mem = reinterpret_cast<intx::uint256*>(&(this->memory->data[ret_offset]));

    //std::cout << "int_mul " << intx::to_string(*a) << " * " << intx::to_string(*b);
	*ret_mem = (*a * *b ) % two_pow_256;
    //std::cout << " = " << intx::to_string(*ret_mem) << std::endl;
}

void BNAPI::div256(uint32_t &a_offset, uint32_t &b_offset, uint32_t &c_offset, uint32_t &ret_offset) {
	intx::uint256* a = reinterpret_cast<intx::uint256*>(&(this->memory->data[a_offset]));
	intx::uint256* b = reinterpret_cast<intx::uint256*>(&(this->memory->data[b_offset]));
	intx::uint256* ret_remainder_mem = reinterpret_cast<intx::uint256*>(&(this->memory->data[ret_offset]));
	intx::uint256* ret_quotient_mem = reinterpret_cast<intx::uint256*>(&(this->memory->data[c_offset]));

	const auto div_res = udivrem(*a, *b);
	*ret_quotient_mem = div_res.quot;
	*ret_remainder_mem = div_res.rem;
}

uint32_t BNAPI::add256(uint32_t &a_offset, uint32_t &b_offset, uint32_t &ret_offset) {
        intx::uint256* a = reinterpret_cast<intx::uint256*>(&(this->memory->data[a_offset]));
        intx::uint256* b = reinterpret_cast<intx::uint256*>(&(this->memory->data[b_offset]));
        intx::uint256* ret_mem = reinterpret_cast<intx::uint256*>(&(this->memory->data[ret_offset]));

        const auto add_res = add_with_carry(*a, *b);

        *ret_mem = std::get<0>(add_res);

        if (tracing) {
            //std::cout << "add256 " << intx::to_string(*a) << " + " << intx::to_string(*b) << " = " << intx::to_string(*ret_mem) << "\n";
        }

        uint32_t carry = 0;
        if (std::get<1>(add_res) > 0) {
          carry = 1;
        }

        return carry;
}

uint32_t BNAPI::sub256(uint32_t &a_offset, uint32_t &b_offset, uint32_t &ret_offset) {
	intx::uint256* a = reinterpret_cast<intx::uint256*>(&(this->memory->data[a_offset]));
	intx::uint256* b = reinterpret_cast<intx::uint256*>(&(this->memory->data[b_offset]));
	intx::uint256* ret_mem = reinterpret_cast<intx::uint256*>(&(this->memory->data[ret_offset]));

	uint32_t carry = 0;
	if (*a < *b) {
	  carry = 1;
	}

	*ret_mem = *a - *b;

    return carry;
}
*/

/*
void add_mod(intx::uint256 *a, intx::uint256 *b, intx::uint256 *modulus, intx::uint256 *ret_mem) {
	intx::uint512 ret_full = intx::uint512{0,*a} + intx::uint512{0,*b};

	if (ret_full >= *modulus) {
	  ret_full -= intx::uint512{0, *modulus};
	}

	*ret_mem = ret_full.lo;
}
*/

/*
void sub_mod(intx::uint256 *a, intx::uint256 *b, intx::uint256 *modulus, intx::uint256 *ret_mem) {
    if (*a < *b) {
      *ret_mem = (*a + *modulus) - *b;
    } else {
      *ret_mem = *a - *b;
    }
}
*/


/*
void mul_mod(uint64_t *a, uint64_t *b, intx::uint256 *modulus, intx::uint256 *inverse, uint64_t *ret) {
    uint64_t *mod = reinterpret_cast<uint64_t *>(modulus);
    uint64_t *inv = reinterpret_cast<uint64_t *>(inverse);
    montgomery_multiplication_256(a, b, mod, inv, ret);
}
*/

/*
void BNAPI::addmod256(uint32_t &a_offset, uint32_t &b_offset, uint32_t &mod_offset, uint32_t &ret_offset) {
  
	intx::uint256* a = reinterpret_cast<intx::uint256*>(&(this->memory->data[a_offset]));
	intx::uint256* b = reinterpret_cast<intx::uint256*>(&(this->memory->data[b_offset]));
  intx::uint256* modulus = reinterpret_cast<intx::uint256*>(&(this->memory->data[mod_offset]));
	intx::uint256* ret_mem = reinterpret_cast<intx::uint256*>(&(this->memory->data[ret_offset]));

    if (tracing) {
        std::cout << "addmod256 " << intx::to_string(*a) << " + " << intx::to_string(*b);
    }

    add_mod(a, b, modulus, ret_mem);

    if (tracing) {
        std::cout <<  " = " << intx::to_string(*ret_mem) << std::endl;
    }
}
*/


void BNAPI::submod256(uint32_t &a_offset, uint32_t &b_offset, uint32_t &mod_offset, uint32_t &ret_offset) {
    printf("BNAPI::submod256\n.");

    /*
    intx::uint256* a = reinterpret_cast<intx::uint256*>(this->memory->data[a_offset]);
    intx::uint256* b = reinterpret_cast<intx::uint256*>(this->memory->data[b_offset]);
    intx::uint256* modulus = reinterpret_cast<intx::uint256*>(this->memory->data[mod_offset]);
    intx::uint256* ret_mem = reinterpret_cast<intx::uint256*>(this->memory->data[ret_offset]);
    */
    /*
    intx::uint256* a = reinterpret_cast<intx::uint256*>(&(this->memory->data[a_offset]));
    intx::uint256* b = reinterpret_cast<intx::uint256*>(&(this->memory->data[b_offset]));
    intx::uint256* modulus = reinterpret_cast<intx::uint256*>(&(this->memory->data[mod_offset]));
    intx::uint256* ret_mem = reinterpret_cast<intx::uint256*>(&(this->memory->data[ret_offset]));
    */


    intx::uint256* a = reinterpret_cast<intx::uint256*>(&(GlobalMem->data[a_offset]));
    intx::uint256* b = reinterpret_cast<intx::uint256*>(&(GlobalMem->data[b_offset]));
    intx::uint256* modulus = reinterpret_cast<intx::uint256*>(&(GlobalMem->data[mod_offset]));
    intx::uint256* ret_mem = reinterpret_cast<intx::uint256*>(&(GlobalMem->data[ret_offset]));


    printf("BNAPI::submod256 casted args.\n.");

    //if (tracing) {
    //    std::cout << "submod256" << intx::to_string(*a) << " - " << intx::to_string(*b);
    //}

    if (*a < *b) {
      *ret_mem = (*a + *modulus) - *b;
    } else {
      *ret_mem = *a - *b;
    }

    printf("BNAPI::submod256 done calcs.\n");

    //sub_mod(a, b, modulus, ret_mem);

    //if (tracing) {
    //    std::cout << " = " << intx::to_string(*ret_mem) << std::endl;
    //}
}


/*
void BNAPI::mulmodmont256(uint32_t &a_offset, uint32_t &b_offset, uint32_t &mod_offset, uint32_t &inv_offset, uint32_t &ret_offset) {
    uint64_t* a = reinterpret_cast<uint64_t*>(&(this->memory->data[a_offset]));
    uint64_t* b = reinterpret_cast<uint64_t*>(&(this->memory->data[b_offset]));
    uint64_t* modulus = reinterpret_cast<uint64_t*>(&(this->memory->data[mod_offset]));
    uint64_t* inv = reinterpret_cast<uint64_t*>(&(this->memory->data[inv_offset]));
    uint64_t* ret = reinterpret_cast<uint64_t*>(&(this->memory->data[ret_offset]));

    if (tracing) {
        intx::uint256 *a_num = reinterpret_cast<intx::uint256 *>(a);
        intx::uint256 *b_num = reinterpret_cast<intx::uint256 *>(b);
        //std::cout << "f1m_mul " << intx::to_string(*a_num) << " * " << intx::to_string(*b_num);
    }

    //mul_mod(a, b, modulus, inv, ret);
    montgomery_multiplication_256(a, b, modulus, inv, ret);

    if (tracing) {
        intx::uint256 *ret_num = reinterpret_cast<intx::uint256 *>(ret);
        //std::cout << " = " << intx::to_string(*ret_num) << std::endl;
    }
}
*/



