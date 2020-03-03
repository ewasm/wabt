

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
      unsigned char prestateRoot[] = {0xb3, 0xc4, 0x18, 0xcb, 0x00, 0xad, 0x7c, 0x90, 0x71, 0x76, 0xbe, 0x86, 0xa5, 0xa2, 0x17, 0x59, 0xb7, 0x4b, 0xd3, 0x82, 0x8e, 0xd6, 0x2a, 0x1e, 0xa2, 0xae, 0x8d, 0xae, 0xa9, 0x8c, 0x5d, 0xa2};

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
      wabt::interp::Memory* mem = env->GetMemory(0);

      unsigned char postStateData[32];

      uint32_t ret_offset = args[0].value.i32;

      uint8_t* mem_ptr = reinterpret_cast<uint8_t*>(&mem->data[ret_offset]);
      uint8_t* mem_ptr_end = reinterpret_cast<uint8_t*>(&mem->data[ret_offset+32]);

      std::copy(mem_ptr, mem_ptr_end, postStateData);


      // print returned state root
      char buffer [65];
      buffer[64] = 0;
      for(int j = 0; j < 32; j++)
        sprintf(&buffer[2*j], "%02X", postStateData[j]);

      std::cout << "eth2_savePostStateRoot: " << std::hex << buffer << std::endl;

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

      wabt::interp::Memory* mem = env->GetMemory(0);

      uint32_t out_offset = args[0].value.i32;
      uint32_t src_offset = args[1].value.i32;
      uint32_t copy_len = args[2].value.i32;

      //std::cout << "eth2_blockDataCopy writing to mem..." << std::endl;
      std::copy(blockData+src_offset, blockData+copy_len, &mem->data[out_offset]);
 
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

      return interp::Result::Ok;
    }
  );


}
