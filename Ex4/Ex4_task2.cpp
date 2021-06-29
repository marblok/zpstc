/*! Laboratory: Advanced signal processing of digital telecommunications
 * (Zaawansowane przetwarzanie sygnałów telekomunikacji cyfrowej)
 *  Ex. 4. task 2.
 *    Sampling rate conversion with fractional delay filters
 *
 * \author Marek Blok
 * \date 2021.06.29
 */
#include <DSP_lib.h>

int main(int argn, char *args[])
{
  /*************************************************************/
  // Log file setup
  DSP::log.SetLogFileName("log_file.txt");
  DSP::log.SetLogState(DSP::e::LogState::file | DSP::e::LogState::console);

  DSP::log << DSP::lib_version_string() << endl << endl;
  /*************************************************************/

/*
 Idea:
 - set of filters is read from a file - L filters designed in MATLAB
 - tests with bandlimited or fullband signals - chirps + spectrograph

 1. Load table of impulse responses of FSD filters

 2. Input -> output buffer with callback

 DSP::u::OutputBuffer (unsigned int BufferSize_in, unsigned int NoOfInputs_in, DSPe_buffer_type cyclic, DSP_clock_ptr ParentClock, DSP_clock_ptr NotificationsClock, unsigned int NoOfOutputs_in, DSP::u::buffer_callback_ptr func_ptr, unsigned int CallbackIdentifier=0)
   If NoOfOutputs_in > 0 NotificationsClock has also the meaning of OutputClock.

 3. Synchronous clocks >> although itmight be better if the asynchronous clocks were used

*/


  DSP::LoadCoef coef_info;
  long int Fp1, Fp2;
  unsigned int L, M;
  unsigned int N_all;
  DSP::Float_vector h_all;

  if (coef_info.Open("ex4_task2_h_FSD_all.coef", "matlab") == false)
  {
  	return -1;
  }
  Fp1 = coef_info.Fp;

  // verify L and calculate M
  DSP::Float_vector tmp;
  coef_info.Load(tmp, 1);
  L = (unsigned int)tmp[0]; M = (unsigned int)tmp[1];
  Fp2 = (Fp1*L)/M;

  // **********************************
  N_all = coef_info.GetSize(0);
  coef_info.Load(h_all, 0);

  /*************************************************************/

  DSP::Clock_ptr InputClock;
  InputClock=DSP::Clock::CreateMasterClock();


  DSP::u::FileInput InputSignal(InputClock, "matlab/test_signal.wav", 1U, DSP::e::SampleType::ST_short, DSP::e::FileType::FT_wav);
  int Fp1_tmp = InputSignal.GetSamplingRate();
  if (Fp1_tmp != Fp1)
  {
    DSP::log << DSP::e::LogMode::Error 
      << "Sampling rate of the input signal is different from sampling rate from fitler coefficients file (ex4_task2_h_FSD_all.coef)" << endl;
  }

  DSP::u::SamplingRateConversion FSDresampler(InputClock, L, M, h_all);

  DSP::u::AudioOutput SoundOut(Fp2, 1, 16);
  DSP::u::FileOutput FileOut_a("ex4_task2.wav", DSP::e::SampleType::ST_short, 1, DSP::e::FileType::FT_wav, Fp2);
  DSP::u::FileOutput FileOut_b("ex4_task2.flt", DSP::e::SampleType::ST_short, 1, DSP::e::FileType::FT_flt, Fp2);

  DSP::log << "Fp1 = " << Fp1 << ", Fp2 = " << Fp2 << ", L = " << L << ", M = " << M << endl;

  /*************************************************************/
  // Connections definitions
  InputSignal.Output("out") >> FSDresampler.Input("in");
  FSDresampler.Output("out") >> SoundOut.Input("in");
  FSDresampler.Output("out") >> FileOut_a.Input("in");
  FSDresampler.Output("out") >> FileOut_b.Input("in");


  /////////////////////////////////
  // check if there are signals
  // connected to all inputs
  DSP::Component::CheckInputsOfAllComponents();

  // *********************************** //
  DSP::Clock::SchemeToDOTfile(InputClock, "Ex4_task2.dot");

  // *********************************** //
  int SamplesInSegment = 512;
  __int64 NoOfSamplesProcessed = 0;
  // 10 seconds
  #define MAX_SAMPLES_TO_PROCESS 10*Fp1
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
          NoOfSamplesProcessed = MAX_SAMPLES_TO_PROCESS - Fp1/5;
    }

    NoOfSamplesProcessed += SamplesInSegment;
    // ********************************************************** //
  }

  /*************************************************************/
  DSP::Clock::FreeClocks();

  return 0;
}
