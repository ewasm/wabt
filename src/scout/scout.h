

#include "src/interp/interp.h"

#include <memory>

#include <iostream>
#include <regex>

#include <fstream>
#include <iomanip>
#include <sstream>
#include <cctype>

using bytes = std::basic_string<uint8_t>;


bytes from_hex(std::string_view hex)
{
    if (hex.length() % 2 == 1)
        //throw std::length_error{"the length of the input is odd"};
        printf("ERROR: the length of the input is odd\n");

    bytes bs;
    int b = 0;
    for (size_t i = 0; i < hex.size(); ++i)
    {
        auto h = hex[i];
        int v;
        if (h >= '0' && h <= '9')
            v = h - '0';
        else if (h >= 'a' && h <= 'f')
            v = h - 'a' + 10;
        else if (h >= 'A' && h <= 'F')
            v = h - 'A' + 10;
        else
            //throw std::out_of_range{"not a hex digit"};
            printf("ERROR: not a hex digit\n");

        if (i % 2 == 0)
            b = v << 4;
        else
            bs.push_back(static_cast<uint8_t>(b | v));
    }
    return bs;
}




using namespace wabt;
using namespace wabt::interp;


void AppendScoutFuncs(wabt::interp::Environment* env, wabt::interp::HostModule* host_module_env) {


  // TODO: read block data from a scout yaml file
  std::ifstream blockDataFile{"./test_block_data.hex"};
  std::string blockdata_hex{std::istreambuf_iterator<char>{blockDataFile}, std::istreambuf_iterator<char>{}};

  blockdata_hex.erase(
      std::remove_if(blockdata_hex.begin(), blockdata_hex.end(), [](char x) { return std::isspace(x); }),
      blockdata_hex.end());


  // bytes is a basic_string
  auto blockdata_bytes = std::make_shared<bytes>(from_hex(blockdata_hex));
  //std::cout << "blockdata bytes length:" << blockdata_bytes->size() << std::endl;

  const unsigned char* blockData = blockdata_bytes->data();

  int block_data_size = std::strlen((char*)blockData);

  //std::cout << "blockData[] length:" << block_data_size << std::endl;
  //printf("blockData: %s\n", blockData);

  //for(int i=0; i < block_data_size; ++i)
  //  std::cout << std::hex << (int)blockData[i];

  //std::cout << std::endl;


  host_module_env->AppendFuncExport(
    "debug_printMemHex",
    {{Type::I32, Type::I32}, {}},
    [env](
      const interp::HostFunc*,
      const interp::FuncSignature*,
      const interp::TypedValues& args,
      interp::TypedValues& results
    ) {
      printf("debug_printMemHex mem_pos: %llu\n", args[0].value.i32);
      return interp::Result::Ok;
    }
  );


  host_module_env->AppendFuncExport(
    "eth2_loadPreStateRoot",
    {{Type::I32}, {}},
    [env](
      const interp::HostFunc*,
      const interp::FuncSignature*,
      const interp::TypedValues& args,
      interp::TypedValues& results
    ) {
      // TODO: read prestateRoot from a scout yaml file

      /***** this is the prestate for the biturbo test case turbo-token-realistic.yaml *****/
      // d2aa91ce2a8031c942f08f7ce4c73843cd6f6a3086649c66ac9bb0528afce5b3
      unsigned char prestateRoot[] = {0xd2, 0xaa, 0x91, 0xce, 0x2a, 0x80, 0x31, 0xc9, 0x42, 0xf0, 0x8f, 0x7c, 0xe4, 0xc7, 0x38, 0x43, 0xcd, 0x6f, 0x6a, 0x30, 0x86, 0x64, 0x9c, 0x66, 0xac, 0x9b, 0xb0, 0x52, 0x8a, 0xfc, 0xe5, 0xb3};

      // printf("eth2_loadPreStateRoot mem_pos: %llu\n", args[0].value.i32);

      uint32_t out_offset = args[0].value.i32;

      wabt::interp::Memory* mem = env->GetMemory(0);

      std::copy(prestateRoot, prestateRoot+32, &mem->data[out_offset]);

      return interp::Result::Ok;
    }
  );


  host_module_env->AppendFuncExport(
    "eth2_savePostStateRoot",
    {{Type::I32}, {}},
    [env](
      const interp::HostFunc*,
      const interp::FuncSignature*,
      const interp::TypedValues& args,
      interp::TypedValues& results
    ) {
      

      //printf("eth2_savePostStateRoot mem_pos: %llu\n", args[0].value.i32);

      wabt::interp::Memory* mem = env->GetMemory(0);

      //printf("eth2_savePostStateRoot got memory.\n");

      unsigned char postStateData[32];

      uint32_t ret_offset = args[0].value.i32;

      uint8_t* mem_ptr = reinterpret_cast<uint8_t*>(&mem->data[ret_offset]);
      uint8_t* mem_ptr_end = reinterpret_cast<uint8_t*>(&mem->data[ret_offset+32]);

      //printf("eth2_savePostStateRoot copying memory...\n");

      std::copy(mem_ptr, mem_ptr_end, postStateData);

      /// print returned state root
      char buffer [33];
      buffer[32] = 0;
      for(int j = 0; j < 16; j++)
        sprintf(&buffer[2*j], "%02X", postStateData[j]);

      std::cout << "eth2_savePostStateRoot: " << std::hex << buffer << std::endl;
      /**** turbo-token test case should print ******
      *     eth2_savePostStateRoot: 3BD1887043B7B8C14B91CCB713A9C900
      */


      return interp::Result::Ok;
    }
  );


  host_module_env->AppendFuncExport(
    "eth2_blockDataSize",
    {{}, {Type::I32}},
    [blockdata_bytes](
      const interp::HostFunc*,
      const interp::FuncSignature*,
      const interp::TypedValues& args,
      interp::TypedValues& results
    ) {
      //printf("eth2_blockDataSize\n");

      int data_size = blockdata_bytes->size();

      //printf("data_size: %d\n", data_size);

      results[0].set_i32(data_size);

      return interp::Result::Ok;
    }
  );



  host_module_env->AppendFuncExport(
    "eth2_blockDataCopy",
    {{Type::I32, Type::I32, Type::I32}, {}},
    [env, blockData](
      const interp::HostFunc*,
      const interp::FuncSignature*,
      const interp::TypedValues& args,
      interp::TypedValues&
    ) {
      //printf("eth2_blockDataCopy.\n");

      // eth2_blockDataCopy(outOffset, srcOffset, length) {
      //wabt::interp::Memory* mem = &env->memories_[0];
      wabt::interp::Memory* mem = env->GetMemory(0);

      //&(mem->data[a_offset])

      // void memorySet(size_t offset, uint8_t value) override { m_wasmMemory->data[offset] = static_cast<char>(value); }

      uint32_t out_offset = args[0].value.i32;
      uint32_t src_offset = args[1].value.i32;
      uint32_t copy_len = args[2].value.i32;

      // TODO: out_offset is incrementing by 266 on every call, which is very weird.
      // should be the same on every call (it is on Scout)

      //printf("eth2_blockDataCopy out_offset: %d\n", args[0].value.i32);
      //printf("eth2_blockDataCopy src_offset: %d\n", args[1].value.i32);
      //printf("eth2_blockDataCopy copy_len: %d\n", args[2].value.i32);


      /*
      //&(mem->data[a_offset]);
      //mem->data[out_offset] = static_cast<char>(value);
      unsigned char dataToCopy[copy_len];
      // std::copy(array+1, array+4, b);
      std::copy(blockData+src_offset, blockData+copy_len, dataToCopy);


      // check that dataToCopy is correct...
      char buffer [33];
      buffer[32] = 0;
      for(int j = 0; j < 16; j++)
        sprintf(&buffer[2*j], "%02X", dataToCopy[j]);

      //std::cout << "eth2_blockDataCopy writing this to mem: " << std::hex << buffer << std::endl;
      */


      //std::cout << "eth2_blockDataCopy writing to mem..." << std::endl;
      std::copy(blockData+src_offset, blockData+copy_len, &mem->data[out_offset]);
      //std::cout << "eth2_blockDataCopy wrote to mem." << std::endl;

      /*
      // inspect written memory
      unsigned char writtenToMem[32];
      uint8_t* mem_ptr = reinterpret_cast<uint8_t*>(&mem->data[out_offset]);
      uint8_t* mem_ptr_end = reinterpret_cast<uint8_t*>(&mem->data[out_offset+32]);

      std::copy(mem_ptr, mem_ptr_end, writtenToMem);

      char bufferWrittenMem [33];
      bufferWrittenMem[32] = 0;
      for(int j = 0; j < 16; j++)
        sprintf(&bufferWrittenMem[2*j], "%02X", writtenToMem[j]);

      std::cout << "eth2_blockDataCopy memory after writing:" << std::hex << bufferWrittenMem << std::endl;
      */


      //env->SetBignumStack(args[0].value.i32);
      return interp::Result::Ok;
    }
  );




}

//void set_i32(uint32_t x) { value.i32 = x; }
//void set_i64(uint64_t x) { value.i64 = x; }