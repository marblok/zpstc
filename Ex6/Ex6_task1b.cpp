/*! Laboratory: Advanced signal processing of digital telecommunications
 * (Zaawansowane przetwarzanie sygnałów telekomunikacji cyfrowej)
 *  Ex. 6. task 1. variant b
 *
 * \author Marek Blok
 * \date 2021.09.15
 */
#include <DSP_lib.h>
#include <memory>

DSP::Float_ptr read_buffer = NULL;
int buffer_size, read_buffer_pointer = -1;
int No_of_samples=0;


int main(int argn, char *args[])
{
  UNUSED_ARGUMENT(argn);
  UNUSED_ARGUMENT(args);

    /*! Idea:
    * - three channels - different modulations - the same shaping filter
    *   DPSK, pi/4-QPSK, (MSK), OQPSK, QAM-16, OOK
    * - FFT - 32 ==> direct generation on intermediate frequency (oversampling 32)
    * - available channels 1..15 (for FMT adjacent channels should not be used)
    *   using: 8, 10, 13
    * other channels are set to zero: SetConstInput
    */
  /*************************************************************/
  // Log file setup
  DSP::log.SetLogFileName("log_file.txt");
  DSP::log.SetLogState(DSP::e::LogState::file | DSP::e::LogState::console);

  DSP::log << DSP::lib_version_string() << std::endl << std::endl;
  /*************************************************************/

  long int Fp2, F_symb;
  DSP::LoadCoef coef_info;
  int N_rc;
  unsigned int L;
  DSP::Float_vector h_rc;

  coef_info.Open("ex6_task1_h_rc.coef", "matlab");
  N_rc = coef_info.GetSize(0);
  if (N_rc < 1)
  {
    DSP::log << DSP::e::LogMode::Error << "No ex6_task1_h_rc.coef: aboarding" << std::endl;
    return -1;
  }
  else
  {
    coef_info.Load(h_rc);
    F_symb = coef_info.Fp;
  }
  /*************************************************************/

  DSP::Clock_ptr BitClock, SymbolClock, OutputClock;
  SymbolClock = DSP::Clock::CreateMasterClock();


  //DSP::u::WaveInput AudioIn(MasterClock, "test.wav", ".");
  //F_symb = AudioIn.GetSamplingRate();

  // DFT length == numer of subchannels
  //  for K is and integer power of two the FFT radix-2 can be used
  int K = 32;
  //int K = 16;
  L = K;

  // symbol rate
  F_symb = 1500; // for L = 32 => Fp2 = 48000


  // output sample rate
  Fp2 = K*F_symb;

  DSP::log << "Fsymb = " << int(F_symb) << ", Fp2 = " << int(Fp2) << ", L = " << L << std::endl;

  BitClock    = DSP::Clock::GetClock(SymbolClock, 2, 1);
  OutputClock = DSP::Clock::GetClock(SymbolClock, K, 1);

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
  // First data channel
  DSP::u::FileInput BinData1(BitClock, "Ex6_task1.cpp", 1U, DSP::e::SampleType::ST_bit, DSP::e::FileType::FT_raw);
  BinData1.SetName("Ex6_task1.cpp");
  //DSP::u::PSKencoder PSKencoder1(DSP::e::PSK_type::QPSK_A);

  DSP::u::FileOutput BinData_out1("ex6_task1_bin_a_new.flt", DSP::e::SampleType::ST_float, 1U, DSP::e::FileType::FT_flt, F_symb);
  BinData_out1.SetName("ex6_task1_bin_a_new.flt");
  BinData1.Output("out") >> BinData_out1.Input("in");

  DSP::u::Serial2Parallel SPconv1(BitClock, 2);
  DSP::u::SymbolMapper PSKencoder1(DSP::e::ModulationType::PSK, 2, 0.0); // QPSK_A
  

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
  // second data channel
  DSP::u::FileInput BinData2(BitClock, "Ex6_task2.cpp", 1U, DSP::e::SampleType::ST_bit, DSP::e::FileType::FT_raw);
  BinData2.SetName("Ex6_task2.cpp");
  DSP::u::Serial2Parallel SPconv2(BitClock, 2);
  DSP::u::SymbolMapper PSKencoder2(DSP::e::ModulationType::PSK, 2, DSP::M_PIf/4); // QPSK_B

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
  // third data channel
  DSP::u::FileInput BinData3(BitClock, "Ex6_task1b.cpp", 1U, DSP::e::SampleType::ST_bit, DSP::e::FileType::FT_raw);
  BinData3.SetName("Ex6_task1.cpp");
  DSP::u::Serial2Parallel SPconv3(BitClock, 2);
  DSP::u::SymbolMapper PSKencoder3(DSP::e::ModulationType::PSK, 2, 0.0); // QPSK_A

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
  BinData1.Output("out") >> SPconv1.Input("in");
  SPconv1.Output("out")  >> PSKencoder1.Input("in");

  BinData2.Output("out") >> SPconv2.Input("in");
  SPconv2.Output("out")  >> PSKencoder2.Input("in");

  BinData3.Output("out") >> SPconv3.Input("in");
  SPconv3.Output("out")  >> PSKencoder3.Input("in");
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
  // FFT with polyphase filters (only for used channels)

  int channel1, channel2, channel3;
  // subchannels numbering convention 0, 1, ..., K-1
  channel1 = 1; channel2 = 4; channel3 = 6;

  auto fft = make_shared<DSP::u::FFT>(K);
  for (int ind = 0; ind < K; ind++)
  {
      if ((ind != channel1) && (ind != channel2) && (ind != channel3))
      {
         string text = "in" + to_string(ind+1);
         fft->SetConstInput(text, 0.0, 0.0);
      }
  }

  // connect narrowband single tone generators (ST)
  // channel no 8
  std::string tekst;

  tekst = "in" + to_string(channel1+1);
  PSKencoder1.Output("out") >> fft->Input(tekst);
  DSP::u::FileOutput SymbData1("ex6_task1a_new.flt", DSP::e::SampleType::ST_float, 2U, DSP::e::FileType::FT_flt, F_symb);
  SymbData1.SetName("ex6_task1a_new.flt");
  PSKencoder1.Output("out") >> SymbData1.Input("in");


  // channel no 10
  tekst = "in" + to_string(channel2+1);
  PSKencoder2.Output("out") >> fft->Input(tekst);
  // channel no 13
  tekst = "in" + to_string(channel3+1);
  PSKencoder3.Output("out") >> fft->Input(tekst);

  // polyphase filters
  vector<shared_ptr<DSP::u::FIR> > H_g(K, nullptr);


  // *************************************************** //
  DSP::u::Parallel2Serial OutputBuffer(SymbolClock, K, 2, true);

  for (int ind = 0; ind < K; ind++)
  {
    // Polyphase filters (joint filters' creation and polyphase decomposition of impulse response h_rc):
    //  Filter with impulse response composed of every K-th sample
    //  of impulse response h_rc starting from (K-1)-th sample.
    H_g[ind] = make_shared<DSP::u::FIR>(true, h_rc, (K-1)-ind, K);

    tekst = "out" + to_string(ind+1);
    fft->Output(tekst) >> H_g[ind]->Input("in");

//    tekst = "in" + to_string(2*ind+1);
//    H_g[ind]->Output("out.re") >> OutputBuffer.Input(tekst);
//    tekst = "in" + to_string(2*ind+2);
//    H_g[ind]->Output("out.im") >> OutputBuffer.Input(tekst);
    tekst = "in" + to_string(ind+1);
    H_g[ind]->Output("out") >> OutputBuffer.Input(tekst);
  }


  // Output to the soundcard
  DSP::u::AudioOutput SoundOut(Fp2, 1, 16);
  // Output to the mono 16bit *.wav file
  DSP::u::FileOutput FileOut2a("ex6_task1b.wav", DSP::e::SampleType::ST_short, 2, DSP::e::FileType::FT_wav, Fp2);
  FileOut2a.SetName("ex6_task1b.wav");
  DSP::u::FileOutput FileOut2b("ex6_task1b.flt", DSP::e::SampleType::ST_float, 2, DSP::e::FileType::FT_flt, Fp2);
  FileOut2b.SetName("ex6_task1b.flt");

  OutputBuffer.Output("out1") >> SoundOut.Input("in");
  OutputBuffer.Output("out") >> FileOut2a.Input("in");
  OutputBuffer.Output("out") >> FileOut2b.Input("in");

  /////////////////////////////////
  // check if there are signals
  // connected to all inputs
  DSP::Component::CheckInputsOfAllComponents();

  // *********************************** //
  DSP::Clock::SchemeToDOTfile(SymbolClock, "ex6_task1b.dot");

  // *********************************** //
  int SamplesInSegment = 512;
  __int64 NoOfSamplesProcessed = 0;
  // 10 seconds
  #define MAX_SAMPLES_TO_PROCESS 5*Fp2
  while(NoOfSamplesProcessed < MAX_SAMPLES_TO_PROCESS)
  {

    // ********************************************************** //
    DSP::Clock::Execute(OutputClock, SamplesInSegment);
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
  H_g.clear();
  fft.reset();

  /*************************************************************/
  DSP::Clock::FreeClocks();

  return 0;
}
