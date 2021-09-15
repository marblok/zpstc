#include <sstream>
using namespace std;

#include <DSP_lib.h>
#include <memory>


void test_symbol_sampler()
{
  map<string,shared_ptr<DSP::Component> > blocks;

  DSP::Clock_ptr BitClock, SymbolClock;
  SymbolClock=DSP::Clock::CreateMasterClock();
  long F_symb = 1500;
  int bits_per_symbol = 2;
  BitClock    = DSP::Clock::GetClock(SymbolClock, bits_per_symbol, 1);

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
  DSP::log << "test_symbol_sampler()" << std::endl;
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
  // First data channel
  blocks["BinData1"] = make_shared<DSP::u::FileInput>(BitClock, "ex6_task1.cpp", 1U, DSP::e::SampleType::ST_bit, DSP::e::FileType::FT_raw);
  blocks["BinData1"]->SetName("ex6_task1.cpp");  //DSP::u::Const BinData1(BitClock, DSP::Float(-1.0));
  blocks["BinData_in"] = make_shared<DSP::u::FileOutput>("bits_in.txt", DSP::e::SampleType::ST_bit_text, 1U, DSP::e::FileType::FT_raw);
  blocks["BinData_in"]->SetName("bits_in.txt");
  blocks["BinData1"]->Output("out") >> blocks["BinData_in"]->Input("in");

  //DSP::u::PSKencoder PSKencoder1(DSP::e::PSK_type::QPSK_A);
  blocks["SPconv"] = make_shared<DSP::u::Serial2Parallel>(BitClock, bits_per_symbol);
  blocks["BinData1"]->Output("out") >> blocks["SPconv"]->Input("in");
  blocks["PSKencoder"] = make_shared<DSP::u::SymbolMapper>(DSP::e::ModulationType::PSK, bits_per_symbol, 0.0); // QPSK_A
  blocks["SPconv"]->Output("out") >> blocks["PSKencoder"]->Input("in");

  blocks["SymbData"] = make_shared<DSP::u::FileOutput>("symbols.flt", DSP::e::SampleType::ST_float, 2U, DSP::e::FileType::FT_flt, F_symb);
  blocks["SymbData"]->SetName("symbols.flt");
  blocks["PSKencoder"]->Output("out") >> blocks["SymbData"]->Input("in");

  blocks["PSKdecoder"] = make_shared<DSP::u::SymbolDemapper>(DSP::e::ModulationType::PSK, bits_per_symbol, 0.0); // QPSK_A
  blocks["PSKencoder"]->Output("out") >> blocks["PSKdecoder"]->Input("in");
  blocks["PSconv"] = make_shared<DSP::u::Parallel2Serial>(SymbolClock, bits_per_symbol);
  blocks["PSKdecoder"]->Output("out") >> blocks["PSconv"]->Input("in");

  // SymbolMapper first sends "empty symbol" which is equivalent to bits_per_symbol zero bits
  blocks["delay"] = shared_ptr<DSP::u::Delay>(new DSP::u::Delay(8-bits_per_symbol)); // byte boundary alignment
  blocks["PSconv"]->Output("out")  >> blocks["delay"]->Input("in");

  blocks["BinData_raw"] = make_shared<DSP::u::FileOutput>("decoded_raw.txt", DSP::e::SampleType::ST_bit, bits_per_symbol, DSP::e::FileType::FT_raw);
  blocks["BinData_raw"]->SetName("decoded_raw.txt");
  blocks["PSKdecoder"]->Output("out") >> blocks["BinData_raw"]->Input("in");

  blocks["BinData_decoded"] = make_shared<DSP::u::FileOutput>("decoded.txt", DSP::e::SampleType::ST_bit, 1U, DSP::e::FileType::FT_raw);
  blocks["BinData_decoded"]->SetName("decoded.txt");
  blocks["delay"]->Output("out") >> blocks["BinData_decoded"]->Input("in");
  (dynamic_cast<DSP::u::FileOutput *>(blocks["BinData_decoded"].get()))->SetSkip(8); // skip first 8 additional zero bits

  blocks["BinData_out"] = make_shared<DSP::u::FileOutput>("bits_out.txt", DSP::e::SampleType::ST_bit_text, 1U, DSP::e::FileType::FT_raw);
  blocks["BinData_out"]->SetName("bits_out.txt");
  blocks["PSconv"]->Output("out") >> blocks["BinData_out"]->Input("in");
  (dynamic_cast<DSP::u::FileOutput *>(blocks["BinData_out"].get()))->SetSkip(0);

  blocks["BinData3"] = make_shared<DSP::u::FileOutput>("bits3.flt", DSP::e::SampleType::ST_float, 1U, DSP::e::FileType::FT_raw);
  blocks["BinData3"]->SetName("bits3.flt");
  blocks["PSconv"]->Output("out") >> blocks["BinData3"]->Input("in");
  (dynamic_cast<DSP::u::FileOutput *>(blocks["BinData3"].get()))->SetSkip(0);

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //



  // *********************************** //
  DSP::Clock::SchemeToDOTfile(SymbolClock, "SymbolMapper-test.dot");
  /////////////////////////////////
  // check if there are signals
  // connected to all inputs
  DSP::Component::CheckInputsOfAllComponents();

  // *********************************** //
  int SamplesInSegment = 512;
  __int64 NoOfSamplesProcessed = 0;
  // 10 seconds
  #define MAX_SAMPLES_TO_PROCESS 50000 //1*Fp2
  while(NoOfSamplesProcessed < MAX_SAMPLES_TO_PROCESS)
  {

    // ********************************************************** //
    DSP::Clock::Execute(SymbolClock, SamplesInSegment);
    // ********************************************************** //

//    if (BinData1.GetBytesRead() > 0)
//    {
//        NoOfSamplesProcessed = 0; // Play the whole file
//    }
//    else // Play 200ms more
//    {
//      if (NoOfSamplesProcessed < MAX_SAMPLES_TO_PROCESS - Fp2/5)
//          NoOfSamplesProcessed = MAX_SAMPLES_TO_PROCESS - Fp2/5;
//    }

    NoOfSamplesProcessed += SamplesInSegment;
    // ********************************************************** //

    //sprintf(tekst, "NoOfSamplesProcessed = %i", int(NoOfSamplesProcessed));
    //DSPf_InfoMessage(tekst);
 }

  // *************************************************************
  blocks.clear();

  DSP::Clock::ListComponents();
  // ************************************************************
  DSP::Clock::FreeClocks();

}

void test_parallel_2_serial()
{
  map<string,shared_ptr<DSP::Component> > blocks;

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
  DSP::log << "test_parallel_2_serial()" << std::endl;
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

  DSP::Clock_ptr BitClock, SymbolClock;
  SymbolClock=DSP::Clock::CreateMasterClock();

  unsigned int no_of_symbol_components = 2; //complex signals
  int bits_per_symbol = 3;
  BitClock    = DSP::Clock::GetClock(SymbolClock, bits_per_symbol, 1);

  // First data channel
  blocks["BinData1"] = make_shared<DSP::u::FileInput>(BitClock, "ex6_task1.cpp", no_of_symbol_components, DSP::e::SampleType::ST_bit, DSP::e::FileType::FT_raw);
  blocks["BinData1"]->SetName("ex6_task1.cpp");
  //DSP::u::Const BinData1(BitClock, DSP::Float(-1.0));
  blocks["BinData1_txt"] = make_shared<DSP::u::FileOutput>("bits_in.txt", DSP::e::SampleType::ST_bit_text, no_of_symbol_components, DSP::e::FileType::FT_raw);
  blocks["BinData1"]->Output("out") >> blocks["BinData1_txt"]->Input("in");

  vector<DSP::Float> SP_init;
  SP_init.resize(bits_per_symbol*no_of_symbol_components, 1.0);
  blocks["SPconv"] = make_shared<DSP::u::Serial2Parallel>(BitClock, bits_per_symbol, no_of_symbol_components, SP_init);
  blocks["BinData1"]->Output("out") >> blocks["SPconv"]->Input("in");
  blocks["BinData_SP_txt"] = make_shared<DSP::u::FileOutput>("bits_SP_out.txt", DSP::e::SampleType::ST_bit_text, bits_per_symbol*no_of_symbol_components, DSP::e::FileType::FT_raw);
  blocks["SPconv"]->Output("out") >> blocks["BinData_SP_txt"]->Input("in");

//  blocks["SymbData"] = make_shared<DSP::u::FileOutput>("symbols.flt", DSP::e::SampleType::ST_float, 2U, DSP::e::FileType::FT_flt);
//  PSKencoder1.Output("out") >> blocks["SymbData"]->Input("in");
//
//  blocks["PSKdecoder"] = make_shared<DSP::u::SymbolDemapper>(SymbolClock, DSP::e::ModulationType::PSK, bits_per_symbol, 0.0); // QPSK_A
//  PSKencoder1.Output("out") >> blocks["PSKdecoder"]->Input("in");

  blocks["PSconv"] = make_shared<DSP::u::Parallel2Serial>(SymbolClock, bits_per_symbol, no_of_symbol_components);
  blocks["SPconv"]->Output("out") >> blocks["PSconv"]->Input("in");

  blocks["BinData2"] = make_shared<DSP::u::FileOutput>("bits.txt", DSP::e::SampleType::ST_bit, no_of_symbol_components, DSP::e::FileType::FT_raw);
  blocks["PSconv"]->Output("out") >> blocks["BinData2"]->Input("in");
  (dynamic_cast<DSP::u::FileOutput *>(blocks["BinData2"].get()))->SetSkip(bits_per_symbol); // skip first bits_per_symbol additional zero symbols

  blocks["BinData2_txt"] = make_shared<DSP::u::FileOutput>("bits_out.txt", DSP::e::SampleType::ST_bit_text, no_of_symbol_components, DSP::e::FileType::FT_raw);
  blocks["PSconv"]->Output("out") >> blocks["BinData2_txt"]->Input("in");
  blocks["BinData2b_txt"] = make_shared<DSP::u::FileOutput>("bits_out_b.txt", DSP::e::SampleType::ST_bit_text, no_of_symbol_components, DSP::e::FileType::FT_raw);
  blocks["PSconv"]->Output("out") >> blocks["BinData2b_txt"]->Input("in");
  (dynamic_cast<DSP::u::FileOutput *>(blocks["BinData2b_txt"].get()))->SetSkip(bits_per_symbol); // skip first bits_per_symbol additional zero symbols

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
  DSP::Clock::SchemeToDOTfile(SymbolClock, "Seria2Parallel-test.dot");
  /////////////////////////////////
  // check if there are signals
  // connected to all inputs
  DSP::Component::CheckInputsOfAllComponents();

  // *********************************** //
  int SamplesInSegment = 512;
  __int64 NoOfSamplesProcessed = 0;
  // 10 seconds
  #define MAX_SAMPLES_TO_PROCESS 50000 //1*Fp2
  while(NoOfSamplesProcessed < MAX_SAMPLES_TO_PROCESS)
  {

    // ********************************************************** //
    DSP::Clock::Execute(SymbolClock, SamplesInSegment);
    // ********************************************************** //

    /*
    if (BinData1.GetBytesRead() > 0)
    {
        NoOfSamplesProcessed = 0; // Play the whole file
    }
    else // Play 200ms more
    {
      if (NoOfSamplesProcessed < MAX_SAMPLES_TO_PROCESS - Fp2/5)
          NoOfSamplesProcessed = MAX_SAMPLES_TO_PROCESS - Fp2/5;
    }
    */

    NoOfSamplesProcessed += SamplesInSegment;
    // ********************************************************** //

    //sprintf(tekst, "NoOfSamplesProcessed = %i", int(NoOfSamplesProcessed));
    //DSPf_InfoMessage(tekst);
 }

  /*************************************************************/
  blocks.clear();

  DSP::Clock::ListComponents();
  /*************************************************************/
  DSP::Clock::FreeClocks();

}

int main(int argn, char *args[])
{
  UNUSED_ARGUMENT(argn);
  UNUSED_ARGUMENT(args);

  /*************************************************************/
  // Log file setup
  DSP::log.SetLogFileName("log_file.txt");
  DSP::log.SetLogState(DSP::e::LogState::file | DSP::e::LogState::console);

  DSP::log << DSP::lib_version_string() << std::endl << std::endl;
  /*************************************************************/

  test_symbol_sampler();

  //test_parallel_2_serial();
  return 0;
}
