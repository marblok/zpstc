#include <DSP_lib.h>

int main(void)
{
  DSP::Clock_ptr SymbolClock;
  SymbolClock=DSP::Clock::CreateMasterClock();

  // Tested: DSP_BPSK, DSP_DBPSK, DSP_DEBPSK
  DSP::u::FileInput BinData1(SymbolClock, "Ex6_task1.cpp", 1U, DSP::e::SampleType::ST_bit, DSP::e::FileType::FT_raw);
  DSP::u::PSKencoder PSKencoder1(DSP::e::PSK_type::DEBPSK);
  DSP::u::PSKdecoder PSKdecoder1(DSP::e::PSK_type::DEBPSK);
  DSP::u::Amplifier PhaseRot(DSP::Complex(0.1,-0.8), 1, true);
  DSP::u::FileOutput BinData_out("ex6_task2_text1.dat", DSP::e::SampleType::ST_bit, 1U, DSP::e::FileType::FT_raw);

  DSP::u::FileOutput Data_in("ex6_task2_text_in.flt", DSP::e::SampleType::ST_float, 1U, DSP::e::FileType::FT_flt);
BinData1.Output("out") >> Data_in.Input("in");
  DSP::u::FileOutput Data_en("ex6_task2_text_en.flt", DSP::e::SampleType::ST_float, 2U, DSP::e::FileType::FT_flt);
PhaseRot.Output("out") >> Data_en.Input("in");
  DSP::u::FileOutput Data_de("ex6_task2_text_de.flt", DSP::e::SampleType::ST_float, 1U, DSP::e::FileType::FT_flt);
PSKdecoder1.Output("out") >> Data_de.Input("in");

BinData1.Output("out") >> PSKencoder1.Input("in");
PSKencoder1.Output("out") >> PhaseRot.Input("in");
PhaseRot.Output("out") >> PSKdecoder1.Input("in");
PSKdecoder1.Output("out") >> BinData_out.Input("in");


  // Tested: DSP::e::PSK_type::QPSK_A, DSP::e::PSK_type::QPSK_B
  DSP::u::FileInput BinData2(SymbolClock, "Ex6_task1.cpp", 2U, DSP::e::SampleType::ST_bit, DSP::e::FileType::FT_raw);
  DSP::u::PSKencoder PSKencoder2(DSP::e::PSK_type::QPSK_A);
  DSP::u::PSKdecoder PSKdecoder2(DSP::e::PSK_type::QPSK_A);
  DSP::u::FileOutput BinData_out2("ex6_task2_text2.dat", DSP::e::SampleType::ST_bit, 2U, DSP::e::FileType::FT_raw);

BinData2.Output("out") >> PSKencoder2.Input("in");
PSKencoder2.Output("out") >> PSKdecoder2.Input("in");
PSKdecoder2.Output("out") >> BinData_out2.Input("in");

// *********************************** //
  DSP::Clock::SchemeToDOTfile(SymbolClock, "ex6_task2_test.dot");
    // *********************************** //
  int SamplesInSegment = 512;
  __int64 NoOfSamplesProcessed = 0;
  // 10 seconds
  #define MAX_SAMPLES_TO_PROCESS 1*1000
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

  DSP::Clock::FreeClocks();

  return 0;
}
