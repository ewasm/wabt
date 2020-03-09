/*
 * Copyright 2017 WebAssembly Community Group participants
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <benchmark/benchmark.h>

//#include <intx/intx.hpp>

#include <algorithm>
#include <cassert>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>
#include <chrono>
#include <iostream>

#include "src/binary-reader.h"
#include "src/cast.h"
#include "src/error-formatter.h"
#include "src/feature.h"
#include "src/interp/binary-reader-interp.h"
#include "src/interp/interp.h"
#include "src/literal.h"
#include "src/option-parser.h"
#include "src/resolve-names.h"
#include "src/stream.h"
#include "src/validator.h"
#include "src/wast-lexer.h"
#include "src/wast-parser.h"

#include "src/scout/scout.h"

using namespace wabt;
using namespace wabt::interp;

using chrono_clock = std::chrono::high_resolution_clock;

static int s_verbose;
static const char* s_infile;
static Thread::Options s_thread_options;
static Stream* s_trace_stream;
static bool s_host_print;
static Features s_features;

static std::unique_ptr<FileStream> s_log_stream;
static std::unique_ptr<FileStream> s_stdout_stream;

Environment env;

DefinedModule* module = nullptr;

enum class RunVerbosity {
  Quiet = 0,
  Verbose = 1,
};

static const char s_description[] =
    R"(  read a file in the wasm binary format, and run in it a stack-based
  interpreter.

examples:
  # parse binary file test.wasm, and type-check it
  $ wasm-interp test.wasm

  # parse test.wasm and run all its exported functions
  $ wasm-interp test.wasm --run-all-exports

  # parse test.wasm, run the exported functions and trace the output
  $ wasm-interp test.wasm --run-all-exports --trace

  # parse test.wasm and run all its exported functions, setting the
  # value stack size to 100 elements
  $ wasm-interp test.wasm -V 100 --run-all-exports
)";

static void ParseOptions(int argc, char** argv) {
  OptionParser parser("wasm-interp", s_description);

  parser.AddOption('v', "verbose", "Use multiple times for more info", []() {
    s_verbose++;
    s_log_stream = FileStream::CreateStdout();
  });
  parser.AddHelpOption();
  s_features.AddOptions(&parser);
  parser.AddOption('V', "value-stack-size", "SIZE",
                   "Size in elements of the value stack",
                   [](const std::string& argument) {
                     // TODO(binji): validate.
                     s_thread_options.value_stack_size = atoi(argument.c_str());
                   });
  parser.AddOption('C', "call-stack-size", "SIZE",
                   "Size in elements of the call stack",
                   [](const std::string& argument) {
                     // TODO(binji): validate.
                     s_thread_options.call_stack_size = atoi(argument.c_str());
                   });
  parser.AddOption('t', "trace", "Trace execution",
                   []() { s_trace_stream = s_stdout_stream.get(); });
  parser.AddOption("host-print",
                   "Include an importable function named \"host.print\" for "
                   "printing to stdout",
                   []() { s_host_print = true; });

  parser.AddArgument("filename", OptionParser::ArgumentCount::One,
                     [](const char* argument) { s_infile = argument; });
  parser.Parse(argc, argv);
}



static wabt::Result ReadModule(const char* module_filename,
                               Environment* env,
                               Errors* errors,
                               DefinedModule** out_module) {
  wabt::Result result;
  std::vector<uint8_t> file_data;

  *out_module = nullptr;

  result = ReadFile(module_filename, &file_data);
  if (Succeeded(result)) {
    const bool kReadDebugNames = true;
    const bool kStopOnFirstError = true;
    const bool kFailOnCustomSectionError = true;
    ReadBinaryOptions options(s_features, s_log_stream.get(), kReadDebugNames,
                              kStopOnFirstError, kFailOnCustomSectionError);
    result = ReadBinaryInterp(env, file_data.data(), file_data.size(), options,
                              errors, out_module);

    if (Succeeded(result)) {
      if (s_verbose) {
        env->DisassembleModule(s_stdout_stream.get(), *out_module);
      }
    }
  }
  return result;
}

static interp::Result PrintCallback(const HostFunc* func,
                                    const interp::FuncSignature* sig,
                                    const TypedValues& args,
                                    TypedValues& results) {
  printf("called host ");
  WriteCall(s_stdout_stream.get(), func->module_name, func->field_name, args,
            results, interp::Result::Ok);
  return interp::Result::Ok;
}

// TODO: this is duplicated in benchmark-interp.cc
// function not implemented here, just need a stub so the import is valid
// the implementation is in interp.cc under `case Opcode::EwasmOpcodeName`
static interp::Result EwasmHostFunc(const HostFunc* func,
                                    const interp::FuncSignature* sig,
                                    const TypedValues& args,
                                    TypedValues& results) {
  return interp::Result::Ok;
}

static void InitEnvironment(Environment* env) {
  if (s_host_print) {
    HostModule* host_module = env->AppendHostModule("host");
    host_module->on_unknown_func_export =
        [](Environment* env, HostModule* host_module, string_view name,
           Index sig_index) -> Index {
      if (name != "print") {
        return kInvalidIndex;
      }

      std::pair<HostFunc*, Index> pair =
          host_module->AppendFuncExport(name, sig_index, PrintCallback);
      return pair.second;
    };
  }

  HostModule* host_module_env = env->AppendHostModule("env");

  // this are here only to make the imports validate
  host_module_env->AppendFuncExport("bignum_f1m_mul", {{Type::I32, Type::I32, Type::I32}, {}}, EwasmHostFunc);
  host_module_env->AppendFuncExport("bignum_f1m_square", {{Type::I32, Type::I32}, {}}, EwasmHostFunc);
  host_module_env->AppendFuncExport("bignum_f1m_add", {{Type::I32, Type::I32, Type::I32}, {}}, EwasmHostFunc);
  host_module_env->AppendFuncExport("bignum_f1m_sub", {{Type::I32, Type::I32, Type::I32}, {}}, EwasmHostFunc);
  host_module_env->AppendFuncExport("bignum_f1m_toMontgomery", {{Type::I32, Type::I32}, {}}, EwasmHostFunc);
  host_module_env->AppendFuncExport("bignum_f1m_fromMontgomery", {{Type::I32, Type::I32}, {}}, EwasmHostFunc);
  host_module_env->AppendFuncExport("bignum_int_mul", {{Type::I32, Type::I32, Type::I32}, {}}, EwasmHostFunc);
  //host_module_env->AppendFuncExport("bignum_int_add", {{Type::I32, Type::I32, Type::I32}, {Type::I32}}, EwasmHostFunc);
  //host_module_env->AppendFuncExport("bignum_int_sub", {{Type::I32, Type::I32, Type::I32}, {Type::I32}}, EwasmHostFunc);
  host_module_env->AppendFuncExport("bignum_int_div", {{Type::I32, Type::I32, Type::I32, Type::I32}, {}}, EwasmHostFunc);

  host_module_env->AppendFuncExport("bignum_add256", {{Type::I32, Type::I32, Type::I32}, {}}, EwasmHostFunc);
  host_module_env->AppendFuncExport("bignum_sub256", {{Type::I32, Type::I32, Type::I32}, {}}, EwasmHostFunc);

  // the scout functions aren't implemented using the "parse call to host func as an opcode"
  // optimization that we use for the bignum host functions. so they're handled in scout.h
  AppendScoutFuncs(env, host_module_env);
}




/*
static wabt::Result ReadAndRunModule(const char* module_filename) {
  constexpr auto to_us = [](chrono_clock::duration d) {
		return std::chrono::duration_cast<std::chrono::microseconds>(d).count();
	};

  wabt::Result result;
  Environment env;
  InitEnvironment(&env);

  Errors errors;
  DefinedModule* module = nullptr;

  const auto parseStartTime = chrono_clock::now();

  result = ReadModule(module_filename, &env, &errors, &module);
  FormatErrorsToFile(errors, Location::Type::Binary);

  const auto now = chrono_clock::now();
	const auto parseDuration = now - parseStartTime;
  const auto execStartTime = now;

  if (Succeeded(result)) {
    Executor executor(&env, s_trace_stream, s_thread_options);
    //RunAllExports(module, &executor, RunVerbosity::Verbose);

    TypedValues args;
    interp::Export* main_ = module->GetExport("main");
    ExecResult exec_result = executor.RunExport(main_, args);

    const auto execFinishTime = chrono_clock::now();
    const auto execDuration = execFinishTime - execStartTime;
    std::cout << "parse time: " << to_us(parseDuration) << "us\n";
    std::cout << "exec time: " << to_us(execDuration) << "us\n";
  }
  return result;
}
*/


static wabt::Result InstantiateModule(const char* module_filename) {
	/*
  constexpr auto to_us = [](chrono_clock::duration d) {
		return std::chrono::duration_cast<std::chrono::microseconds>(d).count();
	};
	*/

  wabt::Result result;
  //Environment env;
  InitEnvironment(&env);

  Errors errors;
  //DefinedModule* module = nullptr;

  //const auto parseStartTime = chrono_clock::now();

  result = ReadModule(module_filename, &env, &errors, &module);
  FormatErrorsToFile(errors, Location::Type::Binary);

  Executor executor(&env, s_trace_stream, s_thread_options);
  ExecResult start_result = executor.RunStartFunction(module);
  if (start_result.result != interp::Result::Ok) {
    WriteResult(s_stdout_stream.get(), "error running start function",
                start_result.result);
    return wabt::Result::Error;;
  }
  //const auto now = chrono_clock::now();
	//const auto parseDuration = now - parseStartTime;
  //std::cout << "parse time: " << to_us(parseDuration) << "us\n";

  return result;
}


// this is called by the benchmark function, so cannot take any args (because the benchmark function can't take any args)
static wabt::Result ExecuteModule() {
  Executor executor(&env, s_trace_stream, s_thread_options);
  TypedValues args;
  interp::Export* main_ = module->GetExport("main");
  ExecResult exec_result = executor.RunExport(main_, args);
  if (exec_result.result == interp::Result::Ok) {
    return wabt::Result::Ok;
  }
  return wabt::Result::Error;
}



/*
int ProgramMain() {
  wabt::Result result = InstantiateModule(s_infile);
  if (Succeeded(result)) {
    wabt::Result result = ExecuteModule();
  }
  return result != wabt::Result::Ok;
}
*/



/*
int ProgramMain(int argc, char** argv) {
  InitStdio();
  s_stdout_stream = FileStream::CreateStdout();

  ParseOptions(argc, argv);

  wabt::Result result = ReadAndRunModule(s_infile);
  return result != wabt::Result::Ok;
}

int main(int argc, char** argv) {
  WABT_TRY
  return ProgramMain(argc, argv);
  WABT_CATCH_BAD_ALLOC_AND_EXIT
}
*/

using namespace benchmark;

namespace
{

  // the benchmark function cannot take any inputs except `state`, which is provided by the benchmark lib
  void wabt_interp(State& state) noexcept
  {
    //printf("wabt_interp...\n");
    // TODO: do we need to reset environment/memory before running benchmark?
    // or reset on every loop?
    wabt::Result result;

    // TODO: InstantiateModule before execute

    for (auto _ : state) {
      //result = InstantiateModule(s_infile);
      result = ExecuteModule();
    }

    //return;
    /*
    // report instantiation and execution time separately using counters??
    auto total_gas_used = int64_t{0};
    auto iteration_gas_used = int64_t{0};
    for (auto _ : state)
        total_gas_used += iteration_gas_used = execute(external_code, external_input);

    state.counters["gas_used"] = Counter(iteration_gas_used);
    state.counters["gas_rate"] = Counter(total_gas_used, Counter::kIsRate);
    */
  }

} // namespace benchmark


int main(int argc, char** argv)
{
    constexpr auto to_us = [](chrono_clock::duration d) {
      return std::chrono::duration_cast<std::chrono::microseconds>(d).count();
    };
    // Initialize is a benchmark.h function?
    Initialize(&argc, argv);

    InitStdio();
    s_stdout_stream = FileStream::CreateStdout();
    ParseOptions(argc, argv);

    const auto parseStartTime = chrono_clock::now();
    wabt::Result result = InstantiateModule(s_infile);
    if (Succeeded(result)) {
      const auto now = chrono_clock::now();
      const auto parseDuration = now - parseStartTime;
      const auto execStartTime = now;
      printf("parse succeeded..\n");
      // test execution before running benchmark..
      result = ExecuteModule();
      if (Succeeded(result)) {
        const auto execFinishTime = chrono_clock::now();
        const auto execDuration = execFinishTime - execStartTime;
        std::cout << "parse time: " << std::dec << to_us(parseDuration) << "us\n";
        std::cout << "exec time: " << std::dec << to_us(execDuration) << "us\n";
        printf("execution succeeded...\n");
        // TODO: do we need to reset environment/memory before running benchmark?
        // or reset on every loop?
        printf("register benchmark...\n");
        RegisterBenchmark("wabt_interp", wabt_interp)->Unit(kMicrosecond);
        printf("run benchmark...\n");
        RunSpecifiedBenchmarks();
        return 0;
      }
    }
    return result != wabt::Result::Ok;
}



