/*! Laboratory: Advanced signal processing of telecommunication signals
 * (Zaawansowane przetwarzanie sygnałów telekomunikacji cyfrowej)
 *  Ex. 1. task 3.
 *   - impulse response loaded from file
 *    (data generated with script: ./matlab/design_task_2.m)
 *   - classical interpolator L = 2
 *
 * \author Marek Blok
 * \date 2021.04.09
 */
#include <DSP_lib.h>

int main(int argn, char *args[])
{
  /*************************************************************/
  // Log file setup
  DSP::log.SetLogFileName("log_file.txt");
  DSP::log.SetLogState(DSP::e::LogState::file | DSP::e::LogState::console);

  DSP::log << DSP::lib_version_string() << std::endl << std::endl;
  /*************************************************************/

  DSP::LoadCoef coef_info;
  int N_LPF;
  DSP::Float_vector h_LPF;

  // loading coefficients from file matlab/ex1_task2.coef
  coef_info.Open("ex1_task2.coef", "matlab");
  N_LPF = coef_info.GetSize(0); // read number of coefficients of first vector (index starts form zero)
  if (N_LPF < 1)
  {
    DSP::log << DSP::e::LogMode::Error << "No filter coefficients: aborting" << std::endl;
    return -1;
  }
  else
  {
    // Load coefficients to h_LPF vector
    coef_info.Load(h_LPF);
  }
  /*************************************************************/

  DSP::Clock_ptr MasterClock;
  // Get pointer to algorithm's master clock
  MasterClock=DSP::Clock::CreateMasterClock();

  long int Fp1, Fp2;

  // Source of samples: file ./DSPElib.wav (current folder)
  DSP::u::WaveInput AudioIn(MasterClock, "DSPElib.wav", ".");
  // get sampling rate associated with block AudioIn
  Fp1 = AudioIn.GetSamplingRate();

  Fp2 = 2*Fp1;
  // two-fold zeroinserter with input working at rate of MasterClock
  DSP::u::Zeroinserter Zeroinserter(MasterClock, 2U);
  // example of appending text to the default block's name (for file *.dot)
  Zeroinserter.SetName("x2");

  // FIR filter with real valued impulse response coefficients (here an interpolation filter)
  // uses samples from vector h_LPF
  DSP::u::FIR InterpFIR(h_LPF);

  // output to the soundcard - sampling rate Fp2 (mono/16bit)
  DSP::u::AudioOutput SoundOut(Fp2, 1, 16);
  // write to file *.wav (mono/16bit, sampling rate Fp2)
  DSP::u::FileOutput FileOut("ex1_task3.wav", DSP::e::SampleType::ST_short, 1, DSP::e::FileType::FT_wav, Fp2);

  /*************************************************************/
  // Blocks connection definitions
  // connecting output "out" of block AudioIn with input "in" of Zeroinserter block
  AudioIn.Output("out") >> Zeroinserter.Input("in");
  // connecting output "out" of Zeroinserter to input "in" of InterpFIR
  Zeroinserter.Output("out") >> InterpFIR.Input("in");
  // connecting output "out" of InterpFIR to input "in" of SoundOut
  InterpFIR.Output("out") >> SoundOut.Input("in");
  // connecting output "out" of InterpFIR to input "in" of FileOut
  InterpFIR.Output("out") >> FileOut.Input("in");
  // Note: one output can be connected to many inputs but many outputs cannot be conneted to one input
  

  /////////////////////////////////
  // Checking if all inputs of all blocks have attached input signals
  DSP::Component::CheckInputsOfAllComponents();

  // *********************************** //
  // Saving the scheme of the implemented algorithm to the * .dot file
  // (this file can be converted to *.gif using script rundot.bat
  //  adapt paths in batch file if necessary)
  DSP::Clock::SchemeToDOTfile(MasterClock, "Ex1_task3.dot");

  // *********************************** //
  // processing in batches of SamplesInSegment input samples
  int SamplesInSegment = 512;

  // counter storing number of input samples processed
  int64_t NoOfSamplesProcessed = 0;

  // Maksimum number of input samples to process (200 miliseconds)
  #define MAX_SAMPLES_TO_PROCESS Fp1/5

  while(NoOfSamplesProcessed < MAX_SAMPLES_TO_PROCESS)
  {

    // ********************************************************** //
    // run the implemented algorithm for the SamplesInSegment of the MasterClock clock cycles
    DSP::Clock::Execute(MasterClock, SamplesInSegment);
    // ********************************************************** //

    if (AudioIn.GetBytesRead() > 0)
    { // reset processed samples counter on successful input file read
      // - forces processing of the whole input file 
      NoOfSamplesProcessed = 0;
    }

    NoOfSamplesProcessed += SamplesInSegment;
    // ********************************************************** //
  }

  /*************************************************************/
  // Releasing reserved resources
  DSP::Clock::FreeClocks();
  /*************************************************************/

  return 0;
}
