#include "intx/intx.hpp"
#include "src/interp/interp.h"

using namespace wabt;
using namespace wabt::interp;

class BNAPI {
private:
    wabt::interp::Memory *memory;
    interp::Result bnHostFunc(const HostFunc* func,
                                     const interp::FuncSignature* sig,
                                     const TypedValues& args,
                                     TypedValues& results);
public:
    void SetMemory(wabt::interp::Memory *memory);
    void AddHostFunctions(wabt::interp::Environment* env, interp::HostModule *host_module);

    void f1m_add(uint32_t &a_offset, uint32_t &b_offset, uint32_t &ret_offset);
    void f1m_sub(uint32_t &a_offset, uint32_t &b_offset, uint32_t &ret_offset);
    void f1m_mul(uint32_t &a_offset, uint32_t &b_offset, uint32_t &ret_offset);
    void f1m_square(uint32_t &a_offset, uint32_t &ret_offset);
    void f1m_fromMont(uint32_t &a_offset, uint32_t &ret_offset);
    void f1m_toMont(uint32_t &a_offset, uint32_t &ret_offset);

    void addmod256(uint32_t &a_offset, uint32_t &b_offset, uint32_t &mod_offset, uint32_t &ret_offset);
    void submod256(uint32_t &a_offset, uint32_t &b_offset, uint32_t &mod_offset, uint32_t &ret_offset);
    void mulmodmont256(uint32_t &a_offset, uint32_t &b_offset, uint32_t &mod_offset, uint32_t &inv_offset, uint32_t &ret_offset);

};
