/*! Laboratory: Advanced signal processing of digital telecommunications
 * (Zaawansowane przetwarzanie sygnałów telekomunikacji cyfrowej)
 *  Ex. 5. task 1.
 *    Sampling rate conversion using I-FIR filter
 *
 * \author Marek Blok
 * \date 2021.07.02
 */
#include <sstream>
#include <iomanip>

#include <DSP_lib.h>

using namespace std;

int main(int argn, char *args[])
{
  /*************************************************************/
  // Log file setup
  DSP::log.SetLogFileName("log_file.txt");
  DSP::log.SetLogState(DSP::e::LogState::file | DSP::e::LogState::console);

  DSP::log << DSP::lib_version_string() << std::endl << std::endl;
  /*************************************************************/

/*
 idea:

 1. Load impulse responses of prototype and image-reject filters

 2. implement a shaping filter + image-reject filter


*/


  DSP::LoadCoef coef_info;
  unsigned int L_IFIR;
  unsigned int N_sh, N_ir;
  DSP::Float_vector tmp;
  DSP::Float_vector h_sh, h_ir;

  // **********************************
  if (coef_info.Open("ex5_task1_h_sh.coef", "matlab") == false)
     return 1;
  N_sh = coef_info.GetSize(0);
  coef_info.Load(h_sh, 0);
  coef_info.Load(tmp, 1);
  L_IFIR = (unsigned int)tmp[0];

  // **********************************
  if (coef_info.Open("ex5_task1_h_ir.coef", "matlab") == false)
  {
     return 1;
  }
  N_ir = coef_info.GetSize();
  if (N_ir == 0)
  {
     return 1;
  }
  coef_info.Load(h_ir);

  /*************************************************************/

  DSP::Clock_ptr InputClock;
  InputClock=DSP::Clock::CreateMasterClock();


  DSP::u::FileInput InputSignal(InputClock, "matlab/delta_44100.wav", 1U, DSP::e::SampleType::ST_short, DSP::e::FileType::FT_wav);
  InputSignal.SetName("delta_44100.wav");
  int Fp1 = InputSignal.GetSamplingRate();
  /*if (Fp1_tmp != Fp1)
  {
    DSPf_ErrorMessage("Problem z sygna�em wej�ciowym");
  } */

  DSP::u::FIR H_sh(h_sh, 0, 1, L_IFIR);
  string H_sh_name = "N_sh = " + to_string(N_sh) + ", L_IFIR = " + to_string(L_IFIR);
  H_sh.SetName(H_sh_name);
  DSP::log << H_sh_name << endl;

  DSP::u::FIR H_ir(h_ir);
  string H_ir_name = "N_ir = " + to_string(N_ir);
  H_ir.SetName(H_ir_name);
  DSP::log << H_ir_name << endl;

//  DSP::u::AudioOutput SoundOut(Fp2, 1, 16);
  DSP::u::FileOutput FileOut_a("matlab/ex5_task1.wav", DSP::e::SampleType::ST_short, 1, DSP::e::FileType::FT_wav, Fp1);
  FileOut_a.SetName("ex5_task1.wav");
  DSP::u::FileOutput FileOut_b("matlab/ex5_task1.flt", DSP::e::SampleType::ST_float, 1, DSP::e::FileType::FT_flt, Fp1);
  FileOut_b.SetName("ex5_task1.flt");

  DSP::u::FileOutput FileOut_test("matlab/ex5_task1a.flt", DSP::e::SampleType::ST_float, 1, DSP::e::FileType::FT_flt, Fp1);
  FileOut_test.SetName("ex5_task1a.flt");

  DSP::log << "Fp1 = " << Fp1 << ", L_IFIR = " << L_IFIR << endl;
  
  /*************************************************************/
  // Connections definitions
  InputSignal.Output("out") >> H_sh.Input("in");
  H_sh.Output("out") >> H_ir.Input("in");
  H_ir.Output("out") >> FileOut_a.Input("in");
  H_ir.Output("out") >> FileOut_b.Input("in");

  H_sh.Output("out") >> FileOut_test.Input("in");

  /////////////////////////////////
  // check if there are signals
  // connected to all inputs
  DSP::Component::CheckInputsOfAllComponents();

  // *********************************** //
  rename("ex5_task1.dot", "ex5_task1.dot~");
  DSP::Clock::SchemeToDOTfile(InputClock, "ex5_task1.dot");

  // *********************************** //
  __int64 start_clk_64, end_clk_64, elapsed_clk_64;

  //clock_t  start_clk, elapsed_clk;
  DSP::Float elapsed_time;



  LARGE_INTEGER lpFrequency, lpPerformanceCount;
  __int64 int64_Freq;
  QueryPerformanceFrequency(&lpFrequency);
  int64_Freq = lpFrequency.QuadPart;
  if (int64_Freq == 0) int64_Freq = 1;

  int SamplesInSegment = 512;
  __int64 NoOfSamplesProcessed = 0;

//  global_stop_watch.Start(); //::wxStartTimer();

//  start_clk=clock();

  QueryPerformanceCounter(&lpPerformanceCount);
  start_clk_64 = lpPerformanceCount.QuadPart;

  // 10 seconds
  __int64 MAX_SAMPLES_TO_PROCESS = 10*Fp1;

  while(NoOfSamplesProcessed < MAX_SAMPLES_TO_PROCESS)
  {

    // ********************************************************** //
    DSP::Clock::Execute(InputClock, SamplesInSegment);
    // ********************************************************** //

    if (InputSignal.GetBytesRead() > 0)
    {
        NoOfSamplesProcessed = 0; // Play the whole file
    }
    else // Play 200ms more
    {
      if (NoOfSamplesProcessed < MAX_SAMPLES_TO_PROCESS - Fp1/5)
          MAX_SAMPLES_TO_PROCESS = NoOfSamplesProcessed + Fp1/5;
    }

    NoOfSamplesProcessed += SamplesInSegment;
    // ********************************************************** //
  }

  QueryPerformanceCounter(&lpPerformanceCount); end_clk_64 = lpPerformanceCount.QuadPart;
  elapsed_clk_64 = end_clk_64-start_clk_64;
  elapsed_time=((DSP::Float)elapsed_clk_64)/int64_Freq; // /CLK_TCK;
    // ********************************************************** //

//    elapsed_clk= clock()-start_clk;
//    elapsed_time=((DSP_float)elapsed_clk)/CLK_TCK;

  DSP::log  << endl;
  DSP::log << DSP::e::LogMode::pause <<  "Elapsed time=" << fixed << setprecision(3) << elapsed_time << "[s], "
  		  << setprecision(2) << NoOfSamplesProcessed/elapsed_time/1000 << "[kSa/s]" << endl;

  /*************************************************************/
  DSP::Clock::FreeClocks();

  return 0;
}
