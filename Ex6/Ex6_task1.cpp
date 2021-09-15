/*! Laboratory: Advanced signal processing of digital telecommunications
 * (Zaawansowane przetwarzanie sygnałów telekomunikacji cyfrowej)
 * Ex. 6. task 1.
 *
 * \author Marek Blok
 * \date 2021.09.15
 */
#include <DSP_lib.h>
#include <memory>

int main(int argn, char *args[])
{
  UNUSED_ARGUMENT(argn);
  UNUSED_ARGUMENT(args);

    /*! Idea:
    * - three channels are used - different modulations - the same shaping filter
    *   DPSK, pi/4-QPSK, (MSK), OQPSK, QAM-16, OOK
    * - FFT - 32 ==> directly placing at intermediate frequency (oversampling 32)
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
    DSP::log << DSP::e::LogMode::Error << "No ex6_task1_h_rc.coef: aborting" << std::endl;
    return -1;
  }
  else
  {
    coef_info.Load(h_rc);
    F_symb = coef_info.Fp;
  }
  /*************************************************************/
  int bits_per_symbol = 2;

  DSP::Clock_ptr BitClock, SymbolClock, OutputClock;
  SymbolClock = DSP::Clock::CreateMasterClock();
  BitClock    = DSP::Clock::GetClock(SymbolClock, bits_per_symbol, 1);


  //DSP::u::WaveInput AudioIn(MasterClock, "test.wav", ".");
  //F_symb = AudioIn.GetSamplingRate();


  // szybkość symbolowa
  F_symb = 1500;

  // długość transformaty DFT / liczba podkanałów
  //  jeżeli K jest całkowitą potęgą dwójki możliwe jest użycie FFT radix-2
//  int K = 16;
  int K = 32;
  L = K;

  // wyjściowa szybkosć próbkowania
  Fp2 = K*F_symb;


  string tekst = "Fsymb = " + to_string(int(F_symb)) + ", Fp2 = " + to_string(int(Fp2)) + ", L = " + to_string(L);
  DSP::log << tekst << std::endl;

  OutputClock=DSP::Clock::GetClock(SymbolClock, K, 1);

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
  // Pierwszy kanał danych
  DSP::u::FileInput BinData1(BitClock, "Ex6_task1.cpp", 1U, DSP::e::SampleType::ST_bit, DSP::e::FileType::FT_raw);
  BinData1.SetName("Ex6_task1.cpp (binary)");
  //DSP::u::PSKencoder PSKencoder1(DSP::e::PSK_type::QPSK_A);
  DSP::u::Serial2Parallel SPconv1(BitClock, bits_per_symbol);
  DSP::u::SymbolMapper PSKencoder1(DSP::e::ModulationType::PSK, bits_per_symbol, 0.0); // QPSK_A
  PSKencoder1.SetName("QPSK_A");

  // ====================================================================== //
  DSP::u::FileOutput FileOutCh1("ex6_task1_ch1.flt", DSP::e::SampleType::ST_float, 2, DSP::e::FileType::FT_flt, F_symb);
  FileOutCh1.SetName("ex6_task1_ch1.flt");

  PSKencoder1.Output("out") >> FileOutCh1.Input("in");
  // ====================================================================== //

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
  // drugi kanał danych
  DSP::u::FileInput BinData2(BitClock, "Ex6_task2.cpp", 1U, DSP::e::SampleType::ST_bit, DSP::e::FileType::FT_raw);
  BinData2.SetName("Ex6_task2.cpp (binary)");
  //DSP::u::PSKencoder PSKencoder2(DSP::e::PSK_type::QPSK_B);
  DSP::u::Serial2Parallel SPconv2(BitClock, bits_per_symbol);
  DSP::u::SymbolMapper PSKencoder2(DSP::e::ModulationType::PSK, bits_per_symbol, DSP::M_PIf/4); // QPSK_B
  PSKencoder2.SetName("QPSK_B");

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
  // trzeci kanał danych
  DSP::u::FileInput BinData3(BitClock, "Ex6_task1b.cpp", 1U, DSP::e::SampleType::ST_bit, DSP::e::FileType::FT_raw);
  BinData3.SetName("Ex6_task1b.cpp (binary)");
  //DSP::u::PSKencoder PSKencoder3(DSP::e::PSK_type::QPSK_A);
  DSP::u::Serial2Parallel SPconv3(BitClock, bits_per_symbol);
  DSP::u::SymbolMapper PSKencoder3(DSP::e::ModulationType::PSK, bits_per_symbol, 0.0); // QPSK_A
  PSKencoder3.SetName("QPSK_A");

  //BinData1.Output("out"), PSKencoder1.Input("in"));
  BinData1.Output("out") >> SPconv1.Input("in");
  SPconv1.Output("out") >> PSKencoder1.Input("in");

  //BinData2.Output("out"), PSKencoder2.Input("in"));
  BinData2.Output("out") >> SPconv2.Input("in");
  SPconv2.Output("out") >> PSKencoder2.Input("in");

  //BinData3.Output("out"), PSKencoder3.Input("in"));
  BinData3.Output("out") >> SPconv3.Input("in");
  SPconv3.Output("out") >> PSKencoder3.Input("in");

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
  // FFT z filtrami polifazowymi (tylko dla używanych kanałów)

  int channel1, channel2, channel3;
  // numery aktywnych podkanałów
  // mumeracja podkanałów: 0, 1, ..., K-1
  channel1 = 1; channel2 = 4; channel3 = 6;

  auto fft = make_shared<DSP::u::FFT>(K);
  for (int ind = 0; ind < K; ind++)
  {
      if ((ind != channel1) && (ind != channel2) && (ind != channel3))
      {
         // definiujemy wejścia o stałej wartości 0+j0 (wszystkie poza używanymi)
         string text = "in" + to_string(ind+1);
         fft->SetConstInput(text, 0.0, 0.0);
      }
  }

  // podłącz generatory wąskopasmowe ST
  // kanał nr 8
  tekst = "in" + to_string(channel1+1);
  PSKencoder1.Output("out") >> fft->Input(tekst);
  // kanał nr 10
  tekst = "in" + to_string(channel2+1);
  PSKencoder2.Output("out") >> fft->Input(tekst);
  // kanał nr 13
  tekst = "in" + to_string(channel3+1);
  PSKencoder3.Output("out") >> fft->Input(tekst);

  // filtry polifazowe
  vector<shared_ptr<DSP::u::FIR> > H_g(K, nullptr);
  vector<shared_ptr<DSP::u::Zeroinserter> > ZeroIns(K, nullptr);
  vector<shared_ptr<DSP::u::Delay> > Z1(K, nullptr);
  vector<shared_ptr<DSP::u::Addition> > Add(K, nullptr);
  for (int ind = 0; ind < K; ind++ )
  {
    // Filtry polifazowe (tworzenie filtrów zintegrowane z dekompozycją polifazową odpowiedzi impulsowej h_rc):
    //  Filtr o odpowiedzi impulsowej złożonej z co K-tej próbki
    //  odpowiedzi impulsowej h_rc zaczynając od próbki (K-1)-ind.
    H_g[ind] = make_shared<DSP::u::FIR>(true, h_rc, (K-1)-ind, K);
    ZeroIns[ind] = make_shared<DSP::u::Zeroinserter>(true, SymbolClock, K);

    string text = "out" + to_string(ind+1);
    fft->Output(text) >> H_g[ind]->Input("in");
    H_g[ind]->Output("out") >> ZeroIns[ind]->Input("in");
  }

  for (int ind = 0; ind <K-1; ind++)
  {
    Z1[ind] = make_shared<DSP::u::Delay>(1,2); // zespolony (!!! potem wariant rzeczywisty)
    Z1[ind]->SetName("Z^-1",false);
    Add[ind] = make_shared<DSP::u::Addition>(0,2); // zespolony (!!! potem wariant rzeczywisty)
  }

  ZeroIns[0]->Output("out") >> Z1[0]->Input("in");
  for (int ind = 1; ind < K; ind++)
  {
     ZeroIns[ind]->Output("out") >> Add[ind-1]->Input("in1");
     Z1[ind-1]->Output("out") >> Add[ind-1]->Input("in2");

     if (ind < K-1)
     {
       Add[ind-1]->Output("out") >> Z1[ind]->Input("in");
     }
  }


  // Output to the soundcard
  DSP::u::AudioOutput SoundOut(Fp2, 1, 16);
  // Output to the mono 16bit *.wav file
  DSP::u::FileOutput FileOut2a("ex6_task1.wav", DSP::e::SampleType::ST_short, 2, DSP::e::FileType::FT_wav, Fp2);
  FileOut2a.SetName("ex6_task1.wav");
  DSP::u::FileOutput FileOut2b("ex6_task1.flt", DSP::e::SampleType::ST_float, 2, DSP::e::FileType::FT_flt, Fp2);
  FileOut2b.SetName("ex6_task1.flt");

  Add[K-2]->Output("out.re") >> SoundOut.Input("in");
  Add[K-2]->Output("out") >> FileOut2a.Input("in");
  Add[K-2]->Output("out") >> FileOut2b.Input("in");

  /////////////////////////////////
  // check if there are signals
  // connected to all inputs
  DSP::Component::CheckInputsOfAllComponents();

  // *********************************** //
  DSP::Clock::SchemeToDOTfile(SymbolClock, "Ex6_task1.dot");

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
  }

  /*************************************************************/
  H_g.clear();
  ZeroIns.clear();
  Z1.clear();
  Add.clear();

  fft.reset();

  /*************************************************************/
  DSP::Clock::FreeClocks();

  return 0;
}
