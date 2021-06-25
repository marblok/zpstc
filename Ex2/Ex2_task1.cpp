/*! Laboratory: Advanced signal processing of digital telecommunications
 * (Zaawansowane przetwarzanie sygnałów telekomunikacji cyfrowej)
 *  Ex. 2. task 1.
 *   - reading filter impulse response from file 
 *    (script generating data file: ./matlab/design_task_1.m
 *   - poliphase interpolator L = 2
 *
 * \author Marek Blok 
 * \date 2021.06.25
 */
#include <DSP_lib.h>

void Processing(void)
{
  DSP::LoadCoef coef_info;
  int N_LPF;
  DSP::Float_vector h_LPF;

  coef_info.Open("ex2_task1.coef", "matlab");
  N_LPF = coef_info.GetSize(0);
  if (N_LPF < 1)
  {  
    DSP::log << DSP::e::LogMode::Error << "No filter coeeficients: aborting" << endl;
    return;
  }
  else
  {
    coef_info.Load(h_LPF);
  }
  /*************************************************************/

  DSP::Clock_ptr MasterClock; 
  MasterClock=DSP::Clock::CreateMasterClock();

  long int Fp1, Fp2;

  DSP::u::WaveInput AudioIn(MasterClock, "DSPElib.wav", ".");
  Fp1 = AudioIn.GetSamplingRate();

  int N_0, N_1;
  DSP::Float_vector g0_LPF((N_LPF+1)/2), g1_LPF((N_LPF+1)/2);
  
  N_0 = 0;
  for (int ind = 0; ind < N_LPF; ind+= 2)
  {
      g0_LPF[N_0] = h_LPF[ind];
      N_0++;
  }
  N_1 = 0;
  for (int ind = 1; ind < N_LPF; ind+= 2)
  {
      g1_LPF[N_1] = h_LPF[ind];
      N_1++;
  }
      
  
  DSP::u::FIR InterpFIR_0(g0_LPF);     
  InterpFIR_0.SetName("0");
  DSP::u::FIR InterpFIR_1(g1_LPF);     
  InterpFIR_1.SetName("1");

  Fp2 = 2*Fp1;
  DSP::u::Zeroinserter Zeroinserter_0(MasterClock, 2U);
  Zeroinserter_0.SetName("0");
  DSP::u::Zeroinserter Zeroinserter_1(MasterClock, 2U);
  Zeroinserter_1.SetName("1");
  
  
  DSP::u::Delay Z_1(1);
  Z_1.SetName("1");
  DSP::u::Addition Sum(2);
  
  // Output to the soundcard 
  DSP::u::AudioOutput SoundOut(Fp2, 1, 16);
  // Output to the mono 16bit *.wav file 
  DSP::u::FileOutput FileOut("ex2_task1.wav", DSP::e::SampleType::ST_short, 1, DSP::e::FileType::FT_wav, Fp2);

  /*************************************************************/
  // Connections definitions
  AudioIn.Output("out") >> InterpFIR_0.Input("in");
  InterpFIR_0.Output("out") >> Zeroinserter_0.Input("in");

  AudioIn.Output("out") >> InterpFIR_1.Input("in");
  InterpFIR_1.Output("out") >> Zeroinserter_1.Input("in");
  Zeroinserter_1.Output("out") >> Z_1.Input("in");

  Zeroinserter_0.Output("out") >> Sum.Input("in1");
  Z_1.Output("out") >> Sum.Input("in2");
  
  Sum.Output("out") >> SoundOut.Input("in");
  Sum.Output("out") >> FileOut.Input("in");
  
  
  /////////////////////////////////
  // check if there are signals 
  // connected to all inputs  
  DSP::Component::CheckInputsOfAllComponents();
  
  // *********************************** //
  DSP::Clock::SchemeToDOTfile(MasterClock, "Ex2_task1.dot");
 
  // *********************************** //
  int SamplesInSegment = 512;
  __int64 NoOfSamplesProcessed = 0;
  // 10 seconds
  #define MAX_SAMPLES_TO_PROCESS 10*Fp1 
  while(NoOfSamplesProcessed < MAX_SAMPLES_TO_PROCESS) 
  {

    // ********************************************************** //
    DSP::Clock::Execute(MasterClock, SamplesInSegment);
    // ********************************************************** //
    
    if (AudioIn.GetBytesRead() > 0)
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

}


int main(int argn, char *args[])
{
  /*************************************************************/
  // Log file setup  
  DSP::log.SetLogFileName("log_file.txt");
  DSP::log.SetLogState(DSP::e::LogState::file | DSP::e::LogState::console);

  DSP::log << DSP::lib_version_string() << endl << endl;
  /*************************************************************/
  
  Processing();

  /*************************************************************/
  DSP::Clock::FreeClocks();
  
  return 0;
}
