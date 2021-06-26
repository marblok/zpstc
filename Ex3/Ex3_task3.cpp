/*! Laboratory: Advanced signal processing of digital telecommunications
 * (Zaawansowane przetwarzanie sygnałów telekomunikacji cyfrowej)
 *  Ex. 3. task 4.
 *   - reading filter impulse response from file 
 *    (script generating data file: ./matlab/design_task_3.m)
 *   - two-stage interpolator L = 20
 *
 * \author Marek Blok
 * \date 2021.06.25
 */
#include <sstream>

#include <DSP_lib.h>

int main(int argn, char *args[])
{
  /*************************************************************/
  // Log file setup
  DSP::log.SetLogFileName("log_file.txt");
  DSP::log.SetLogState(DSP::e::LogState::file | DSP::e::LogState::console);

  DSP::log << DSP::lib_version_string() << endl << endl;
  /*************************************************************/

  long int Fp1, Fp2, F_symb;
  DSP::LoadCoef coef_info;
  int N_rc, N2;
  unsigned int L1, L2;
  DSP::Float_vector h_rc, h2;

  coef_info.Open("ex3_task3_h_rc.coef", "matlab");
  N_rc = coef_info.GetSize(0);
  if (N_rc < 1)
  {
    DSP::log << DSP::e::LogMode::Error << "No ex3_task3_h_rc.coef: aborting" << endl;
    return -1;
  }
  else
  {
    coef_info.Load(h_rc);
    Fp1 = coef_info.Fp;
  }
  /*************************************************************/
  coef_info.Open("ex3_task3_h2.coef", "matlab");
  N2 = coef_info.GetSize(0);
  if (N2 < 1)
  {
    DSP::log << DSP::e::LogMode::Error << "No ex3_task3_h2.coef: aborting" << endl;
    return -1;
  }
  else
  {
    coef_info.Load(h2);
    Fp2 = coef_info.Fp;
  }
  /*************************************************************/

  DSP::Clock_ptr SymbolClock, SecondClock;
  SymbolClock=DSP::Clock::CreateMasterClock();


  //DSP::u::WaveInput AudioIn(MasterClock, "test.wav", ".");
  //F_symb = AudioIn.GetSamplingRate();

  DSP::u::FileInput BinData(SymbolClock, "ex3_task3.cpp", 2U, DSP::e::SampleType::ST_bit, DSP::e::FileType::FT_raw);
  F_symb = 2400;
  DSP::u::PSKencoder PSKencoder(DSP::e::PSK_type::QPSK_A);

  L1 = Fp1 / F_symb;
  L2 = Fp2 / Fp1;
  DSP::log << "Fsymb = " << F_symb << ", Fp1 = " << Fp1 << ", Fp2 = " << Fp2 << ", L1 = " << L1 << ", L2 = " << L2 << endl;

  SecondClock=DSP::Clock::GetClock(SymbolClock, L1, 1);

  DSP::u::SamplingRateConversion SRC1(true, SymbolClock, L1, 1, h_rc);
  SRC1.SetName("SRC1");

  DSP::u::SamplingRateConversion SRC2(true, SecondClock, L2, 1, h2);
  SRC2.SetName("SRC2");
  DSP::u::DDScos Heter(SRC2.GetOutputClock(), true, 0.5, (DSP::M_PIx2*2500)/DSP::Float(Fp2));
  DSP::u::Multiplication Mul(0, 2);
  DSP::u::Vacuum V1;


  // Output to the soundcard
  DSP::u::AudioOutput SoundOut(Fp2, 1, 16);
  DSP::u::FileOutput FileOut1("ex3_task3a.flt", DSP::e::SampleType::ST_float, 2, DSP::e::FileType::FT_flt, Fp1);
  // Output to the mono 16bit *.wav file
  DSP::u::FileOutput FileOut2a("ex3_task3b.wav", DSP::e::SampleType::ST_short, 1, DSP::e::FileType::FT_wav, Fp2);
  DSP::u::FileOutput FileOut2b("ex3_task3b.flt", DSP::e::SampleType::ST_float, 1, DSP::e::FileType::FT_flt, Fp2);

  /*************************************************************/
  // Connections definitions
  //AudioIn.Output("out") >> SRC1.Input("in");
  BinData.Output("out") >> PSKencoder.Input("in");
  PSKencoder.Output("out") >> SRC1.Input("in");

  SRC1.Output("out") >> FileOut1.Input("in");
  SRC1.Output("out") >> SRC2.Input("in");

  SRC2.Output("out") >> Mul.Input("in1");
  Heter.Output("out") >> Mul.Input("in2");

  Mul.Output("out.im") >> V1.Input("in");

  Mul.Output("out.re") >> SoundOut.Input("in");
  Mul.Output("out.re") >> FileOut2a.Input("in");
  Mul.Output("out.re") >> FileOut2b.Input("in");


  /////////////////////////////////
  // check if there are signals
  // connected to all inputs
  DSP::Component::CheckInputsOfAllComponents();

  // *********************************** //
  DSP::Clock::SchemeToDOTfile(SymbolClock, "Ex3_task3.dot");

  // *********************************** //
  int SamplesInSegment = 4*512;
  __int64 NoOfSamplesProcessed = 0;
  // 10 seconds
  #define MAX_SAMPLES_TO_PROCESS 10*Fp1
  while(NoOfSamplesProcessed < MAX_SAMPLES_TO_PROCESS)
  {

    // ********************************************************** //
    DSP::Clock::Execute(SymbolClock, SamplesInSegment);
    // ********************************************************** //

    int bytes_read = BinData.GetBytesRead();
    DSP::log << "BinData.GetBytesRead() = " << bytes_read << endl;
    if (bytes_read > 0)
    {
        NoOfSamplesProcessed = 0; // Play the whole file
    }
    else // Play 200ms more
    {
      if (NoOfSamplesProcessed < MAX_SAMPLES_TO_PROCESS - Fp1/5)
          NoOfSamplesProcessed = MAX_SAMPLES_TO_PROCESS - Fp1/5;
    }

    NoOfSamplesProcessed += SamplesInSegment;
    // ********************************************************** //
  }

  /*************************************************************/
  DSP::Clock::FreeClocks();

  return 0;
}
